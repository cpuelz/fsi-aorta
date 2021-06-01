// Filename: VelocityBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <four_chambered_heart/VelocityBcCoefs.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs(const CirculationModel* circ_model, const int comp_idx)
    : d_circ_model(circ_model), d_comp_idx(comp_idx)
{
    // intentionally blank
    return;
} // VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs

void
VelocityBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                            Pointer<ArrayData<NDIM, double> >& bcoef_data,
                            Pointer<ArrayData<NDIM, double> >& gcoef_data,
                            const Pointer<hier::Variable<NDIM> >& /*variable*/,
                            const Patch<NDIM>& patch,
                            const BoundaryBox<NDIM>& bdry_box,
                            double /*fill_time*/) const
{
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif
    const Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);

        // default to zero-pressure BCs:
        a = (axis == d_comp_idx) ? 0.0 : 1.0;
        b = (axis == d_comp_idx) ? 1.0 : 0.0;
        g = 0.0;

        // loop over circulation model sources on this location index
        const int num_sources_here = d_circ_model->d_location_index_to_src_id[location_index].size();
        if(num_sources_here > 0)
        {
            std::vector<double> psrc(num_sources_here);
            std::vector<double> usrc(num_sources_here);
            std::vector<double> rsrc(num_sources_here);
            std::vector<double> r(num_sources_here);
            std::vector<bool> use_velocity_bcs(num_sources_here);
            for (int ss = 0; ss < num_sources_here; ++ss)
            {
                const int source_id = d_circ_model->d_location_index_to_src_id[location_index][ss];
                psrc[ss] = d_circ_model->d_psrc[source_id];
                usrc[ss] = d_circ_model->d_usrc_prescribed[source_id];
                rsrc[ss] = d_circ_model->d_rsrc[source_id];
                use_velocity_bcs[ss] = d_circ_model->d_use_velocity_bcs[source_id];
                IBTK::Point posn;
                for (int d = 0; d < NDIM; ++d)
                {
                    posn(d) = d_circ_model->d_posn[source_id*NDIM + d];
                }
                
                double X[NDIM];
                double r_sq = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                    if (d != axis) r_sq += pow(X[d] - posn[d], 2.0);
                }
                r[ss] = sqrt(r_sq);

            }
            
            // now loop over the sources here and apply the right boundary conditions
            bool in_an_active_source = false;
            static const double delta = 0.2;
            for (int ss = 0; ss < num_sources_here; ++ss)
            {
                if (r[ss] < rsrc[ss])
                {
                    in_an_active_source = true;
                    if(use_velocity_bcs[ss]) // use velocity (flow) boundary conditions
                    {
                        const double umax = 2*usrc[ss];
                        a = 1.0;
                        b = 0.0;
                        if (axis == d_comp_idx)
                        {
                            g = umax * (1.0 - r[ss] * r[ss] / (rsrc[ss] * rsrc[ss])) * ((side == 0) ? 1.0 : -1.0); 
                        }
                        else
                        {
                            g = 0.0;
                        }
                    }
                    else // use pressure boundary conditions
                    {
                        if (axis == d_comp_idx) // impose normal component of fluid traction
                        {
                            a = 0.0;
                            b = 1.0;
                            g = -psrc[ss];
                        }
                        if (axis != d_comp_idx) // set tangential velocities to zero.
                        {
                            a = 1.0;
                            b = 0.0;
                            g = 0.0;
                        }                            
                    }
                }
                else if (r[ss] < rsrc[ss] + delta)
                {
                    a = 1.0;
                    b = 0.0;
                    g = 0.0;
                }
            }
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
