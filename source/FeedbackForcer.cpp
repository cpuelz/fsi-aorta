// Filename: FeedbackForcer.cpp
// Created on 04 May 2007 by Boyce Griffith

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <tbox/Utilities.h>

#include <four_chambered_heart/FeedbackForcer.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
} // namespace

////////////////////////////// PUBLIC ///////////////////////////////////////

FeedbackForcer::FeedbackForcer(const CirculationModel* circ_model,
                               const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
    : d_circ_model(circ_model), d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
{
    // intentionally blank
    return;
} // FeedbackForcer

FeedbackForcer::~FeedbackForcer()
{
    // intentionally blank
    return;
} // ~FeedbackForcer

bool
FeedbackForcer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
FeedbackForcer::setDataOnPatch(const int data_idx,
                               Pointer<Variable<NDIM> > /*var*/,
                               Pointer<Patch<NDIM> > patch,
                               const double /*data_time*/,
                               const bool initial_time,
                               Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(F_data);
#endif
    F_data->fillAll(0.0);
    if (initial_time) return;
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    Pointer<SideData<NDIM, double> > U_current_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    Pointer<SideData<NDIM, double> > U_new_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(U_current_data);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());

    double L[NDIM];
    int offset[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        L[d] = max(dx_coarsest[d], 2.0 * dx_finest[d]);
        offset[d] = static_cast<int>(L[d] / dx[d]);
    }

    // Start with no force anywhere.
    F_data->fillAll(0.0);

    // Clamp all tangential velocities.
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            if (!pgeom->getTouchesRegularBoundary(axis, side)) continue;

            Box<NDIM> bdry_box = domain_box;
            if (side == 0)
            {
                bdry_box.upper(axis) = domain_box.lower(axis) + offset[axis];
            }
            else
            {
                bdry_box.lower(axis) = domain_box.upper(axis) - offset[axis];
            }
            bdry_box = bdry_box * patch_box;

            for (int component = 0; component < NDIM; ++component)
            {
                if (component == axis) continue;
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                    }
                    const double x_bdry = side == 0 ? x_lower[axis] : x_upper[axis];
                    (*F_data)(i_s) = smooth_kernel((X[axis] - x_bdry) / L[axis]) * (-kappa * U);
                }
            }
        }
    }

    // Attempt to prevent flow reversal points near the domain boundary.
    for (int ss = 0; ss < d_circ_model->d_nsrc; ++ss)
    {
        const double qsrc = d_circ_model->d_qsrc[ss];
        const double rsrc = d_circ_model->d_rsrc[ss];
        const int axis = d_circ_model->d_src_axis[ss];
        const int side = d_circ_model->d_src_side[ss];
        const bool is_lower = side == 0;
        IBTK::Point posn;
        for (int d = 0; d < NDIM; ++d)
        {
            posn(d) = d_circ_model->d_posn[ss * NDIM + d];
        }
        const bool inflow_bdry = qsrc < 0.0;
        const bool outflow_bdry = qsrc > 0.0;
        if (pgeom->getTouchesRegularBoundary(axis, side))
        {
            Box<NDIM> bdry_box = domain_box;
            if (side == 0)
            {
                bdry_box.upper(axis) = domain_box.lower(axis) + offset[axis];
            }
            else
            {
                bdry_box.lower(axis) = domain_box.upper(axis) - offset[axis];
            }
            bdry_box = bdry_box * patch_box;
            for (int component = 0; component < NDIM; ++component)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                    }
                    double r_sq = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (d != axis) r_sq += pow(X[d] - posn[d], 2.0);
                    }
                    const double r = sqrt(r_sq);
                    double mask = (component == axis ? 0.0 : 1.0);
                    if (component == axis && r <= rsrc)
                    {
                        const double n = is_lower ? -1.0 : +1.0;
                        const double U_dot_n = U * n;
                        if ((inflow_bdry && U_dot_n > 0.0) || (outflow_bdry && U_dot_n < 0.0))
                        {
                            mask = 1.0;
                        }
                    }
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    mask *= smooth_kernel((X[axis] - x_bdry) / L[axis]);
                    (*F_data)(i_s) = mask * (-kappa * U);
                }
            }
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
