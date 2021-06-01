// Filename: SourceDistributer.cpp
// Created on 02 Jun 2018 by Boyce Griffith

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

// IBAMR INCLUDES
#include <ibtk/IndexUtilities.h>

#include <four_chambered_heart/SourceDistributer.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
cos_kernel(const double x, const double eps)
{
    if (std::abs(x) > eps)
    {
        return 0.0;
    }
    else
    {
        return 0.5 * (1.0 + cos(M_PI * x / eps)) / eps;
    }
} // cos_kernel
} // namespace

////////////////////////////// PUBLIC ///////////////////////////////////////

SourceDistributer::SourceDistributer(std::vector<std::unique_ptr<Source> > sources,
                                     const INSHierarchyIntegrator* fluid_solver,
                                     const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
    : d_sources(std::move(sources)), d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
{
    // intentionally blank
    return;
} // SourceDistributer

SourceDistributer::~SourceDistributer()
{
    // intentionally blank
    return;
} // ~SourceDistributer

bool
SourceDistributer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
SourceDistributer::setDataOnPatch(const int data_idx,
                                  Pointer<hier::Variable<NDIM> > /*var*/,
                                  Pointer<Patch<NDIM> > patch,
                                  const double /*data_time*/,
                                  const bool initial_time,
                                  Pointer<PatchLevel<NDIM> > patch_level)
{
    Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(data_idx);
    TBOX_ASSERT(Q_data);
    Q_data->fillAll(0.0);
    if (initial_time) return;
    const IntVector<NDIM>& ratio = patch_level->getRatio();
    const Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const dx = patch_geom->getDx();
    const double* const x_lower = patch_geom->getXLower();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geom->getDx();
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();

    // The source radius must be an integer multiple of the grid spacing.
    //
    // \todo The source radius should probably be set in the Source data
    // structure.
    boost::array<double, NDIM> r;
    for (int d = 0; d < NDIM; ++d)
    {
        r[d] = 4.0 * dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }

    // Spread the source strength onto the Cartesian grid.
    for (unsigned int n = 0; n < d_sources.size(); ++n)
    {
        // Determine the position of the source.
        const std::vector<double> x_src(&d_sources[n]->d_current_location[0],
                                        &d_sources[n]->d_current_location[0] + NDIM);
        const double q_src = d_sources[n]->d_q_src;

        // Determine the approximate source stencil box.
        const hier::Index<NDIM> i_center = IndexUtilities::getCellIndex(x_src, grid_geom, ratio);
        Box<NDIM> stencil_box(i_center, i_center);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            stencil_box.grow(d, static_cast<int>(r[d] / dx[d]) + 1);
        }

        // Set the Cartesian grid values.
        for (Box<NDIM>::Iterator b(patch_box * stencil_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            double wgt = 1.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double x_center = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                wgt *= cos_kernel(x_center - x_src[d], r[d]);
            }
            (*Q_data)(i) += q_src * wgt;
        }
    }
    return;
}

void
SourceDistributer::interpolatePressure(const int p_data_idx, const double /*data_time*/)
{
    // Compute the mean pressure at the sources/sinks associated with each level
    // of the Cartesian grid.
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geom->getDx();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    // The source radius must be an integer multiple of the grid spacing.
    //
    // \todo The source radius should probably be set in the Source data
    // structure.
    boost::array<double, NDIM> r;
    for (int d = 0; d < NDIM; ++d)
    {
        r[d] = 4.0 * dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }

    // Interpolate the pressure from the Cartesian grid.
    std::vector<double> p_src(d_sources.size(), 0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        Pointer<PatchLevel<NDIM> > finer_level =
            (ln < finest_ln ? d_patch_hierarchy->getPatchLevel(ln + 1) : Pointer<BasePatchLevel<NDIM> >(NULL));
        const IntVector<NDIM>& finer_ratio = (ln < finest_ln ? finer_level->getRatio() : IntVector<NDIM>(1));
        const Box<NDIM> finer_domain_box_level = Box<NDIM>::refine(domain_box, finer_ratio);
        const hier::Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const hier::Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        boost::array<double, NDIM> finer_dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d] / static_cast<double>(finer_ratio(d));
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& patch_lower = patch_box.lower();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            Pointer<CellData<NDIM, double> > P_data = patch->getPatchData(p_data_idx);
            for (unsigned int n = 0; n < d_sources.size(); ++n)
            {
                // Determine the position of the source.
                const std::vector<double>& x_src = d_sources[n]->d_current_location;
                TBOX_ASSERT(x_src.size() == NDIM);

                // Determine the approximate source stencil box.
                const hier::Index<NDIM> i_center = IndexUtilities::getCellIndex(x_src, grid_geom, ratio);
                Box<NDIM> stencil_box(i_center, i_center);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    stencil_box.grow(d, static_cast<int>(r[d] / dx[d]) + 1);
                }

                // Read the Cartesian grid values.
                for (Box<NDIM>::Iterator b(patch_box * stencil_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    double wgt = 1.0;
                    std::vector<double> x_center(NDIM, 0.0);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x_center[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                        wgt *= cos_kernel(x_center[d] - x_src[d], r[d]) * dx[d];
                    }
                    // check to see if center of cell is actually on the finer level.
                    // if it is, then we don't include it in the calculation of the pressure
                    // on this level.
                    const hier::Index<NDIM> finer_i = IndexUtilities::getCellIndex(x_center,
                                                                                   domainXLower,
                                                                                   domainXUpper,
                                                                                   finer_dx.data(),
                                                                                   finer_domain_box_level_lower,
                                                                                   finer_domain_box_level_upper);
                    if (ln == finest_ln || !finer_level->getBoxes().contains(finer_i))
                    {
                        p_src[n] += (*P_data)(i)*wgt;
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&p_src[0], d_sources.size());
    for (unsigned int n = 0; n < d_sources.size(); ++n)
    {
        d_sources[n]->d_p_src = p_src[n];
    }
    return;
}

void
SourceDistributer::updateSources(IBFEMethod* ib_method_ops, const double loop_time, const double dt)
{
    unsigned int source_n = 0;
    for (std::unique_ptr<Source>& source : d_sources)
    {
        source->updateCurrentLocation(ib_method_ops);
        source->updateSourceStrength(loop_time, dt, source_n);
        ++source_n;
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
