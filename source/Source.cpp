/*
 * File:   Source.cpp
 * Author: cpuelz
 *
 * Created on June 2, 2018, 2:48 PM
 */

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/point.h>

#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/Source.h>

namespace
{
    // Conversion factors.
    static const double flconv = 1 / 0.06;  // l/min ===> ml/sec
    static const double prconv = 1333.2239; // mmHg  ===> dyne/cm^2
}

Source::Source(const string& object_name,
               const int part,
               Pointer<Database> source_db,
               IBFEMethod* ib_method_ops,
               bool register_for_restart)
    : d_time(0.0),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_part(part),
      d_initial_location(NDIM, 0.0),
      d_current_location(NDIM, 0.0)
{
#if !defined(NDEBUG)
    assert(!object_name.empty());
#endif
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }
    if (source_db)
    {
        // get source parameters from input file
        d_is_pressure_meter = source_db->getBool("is_pressure_meter");
        d_is_pressure_source = source_db->getBool("is_pressure_source");
        if (!d_is_pressure_source)
        {
            d_prescribed_flow = source_db->getDouble("prescribed_flow");
        }
        d_R_src = source_db->getDouble("R_src");
        d_L_src = source_db->getDouble("L_src");
        d_ramp_time = source_db->getDouble("ramp_time");
        d_P_reservoir = source_db->getDouble("P_reservoir");
        d_sideset_names = source_db->getStringArray("sideset_names");
        d_data_directory_name = source_db->getString("data_directory_name");
        d_data_filename = source_db->getString("data_filename");

        // make directory for storing source data
        Utilities::recursiveMkdir(d_data_directory_name);

        // initialize sideset boundary IDs
        for (unsigned int ii = 0; ii < unsigned(d_sideset_names.size()); ++ii)
        {
            d_sideset_IDs.push_back(MeshInfo::sideset_map[d_part][d_sideset_names[ii]]);
        }

        const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
        const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
        const MeshBase* mesh = &equation_systems->get_mesh();
        const BoundaryInfo& boundary_info = *mesh->boundary_info;

        // compute initial location of fluid source
        const System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
        const unsigned int coords_sys_num = X_system.number();
        NumericVector<double>& X_ghost_vec = *X_system.current_local_solution;
        copy_and_synch(*X_system.solution, X_ghost_vec);

        // loop over elements
        MeshBase::const_element_iterator el = mesh->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = mesh->active_local_elements_end();

        int local_num_elem_sides = 0;
        std::vector<double> centroid(NDIM, 0.0);
        // loop over local elements
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;

            // looping over element sides
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                // get sideset IDs
                std::vector<boundary_id_type> bdry_ids;
                boundary_info.boundary_ids(elem, side, bdry_ids);

                if (find_first_of(bdry_ids.begin(), bdry_ids.end(), d_sideset_IDs.begin(), d_sideset_IDs.end()) !=
                    bdry_ids.end())
                {
                    libMesh::UniquePtr<const Elem> elem_side = elem->build_side_ptr(side);
                    // loop over nodes over this element side
                    for (unsigned int nn = 0; nn < elem_side->n_nodes(); ++nn)
                    {
                        const Node* const node_ptr = elem_side->node_ptr(nn);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            const int coords_dof_idx = node_ptr->dof_number(coords_sys_num, d, 0);
                            centroid[d] += X_ghost_vec.el(coords_dof_idx) / static_cast<double>(elem_side->n_nodes());
                        }
                    }
                    local_num_elem_sides += 1;
                }
            }
        }

        // collect info across all processors
        const int total_num_elem_sides = SAMRAI_MPI::sumReduction(local_num_elem_sides);
        SAMRAI_MPI::sumReduction(&centroid[0], NDIM);

        // finish computing centroid
        for (int d = 0; d < NDIM; ++d)
        {
            centroid[d] /= static_cast<double>(total_num_elem_sides);
            d_initial_location[d] = centroid[d];
        }
    }
    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    else
    {
        // initialize current location
        for (int d = 0; d < NDIM; ++d)
        {
            d_current_location[d] = d_initial_location[d];
        }

        // initialize dynamic variables
        d_q_src = 0.0;
        d_p_src = 0.0;
    }
}

void
Source::updateCurrentLocation(IBFEMethod* ib_method_ops)
{
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;

    const System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    const unsigned int coords_sys_num = X_system.number();
    NumericVector<double>& X_ghost_vec = *X_system.current_local_solution;

    std::vector<double> centroid(NDIM, 0.0);

    // fortunately we only need to compute this once since we never remesh the
    // Lagrangian stuff. Extract all relevant boundary dof ids now and then we
    // can quickly look them up after each time step.
    if (d_ibfe_method != ib_method_ops)
    {
        d_ibfe_method = ib_method_ops;
        for (unsigned int d = 0; d < NDIM; ++d) d_locally_owned_side_dofs[d].resize(0);
        d_local_num_elem_sides = 0;
        d_nodes_per_side = -1;

        MeshBase::const_element_iterator el = mesh->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = mesh->active_local_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;

            // looping over element sides
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                // get sideset IDs
                std::vector<boundary_id_type> bdry_ids;
                boundary_info.boundary_ids(elem, side, bdry_ids);

                if (find_first_of(bdry_ids.begin(), bdry_ids.end(), d_sideset_IDs.begin(), d_sideset_IDs.end()) !=
                    bdry_ids.end())
                {
                    std::unique_ptr<const Elem> elem_side = elem->build_side_ptr(side);
                    // loop over nodes over this element side
                    for (unsigned int nn = 0; nn < elem_side->n_nodes(); ++nn)
                    {
                        const Node* const node_ptr = elem_side->node_ptr(nn);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            const int coords_dof_idx = node_ptr->dof_number(coords_sys_num, d, 0);
                            d_locally_owned_side_dofs[d].push_back(coords_dof_idx);
                        }
                    }
                    if (d_nodes_per_side == -1)
                    {
                        d_nodes_per_side = elem_side->n_nodes();
                    }
                    else
                    {
                        TBOX_ASSERT(d_nodes_per_side == static_cast<int>(elem_side->n_nodes()));
                    }
                    d_local_num_elem_sides += 1;
                }
            }
        }
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        for (const auto coords_dof_idx : d_locally_owned_side_dofs[d])
            centroid[d] += X_ghost_vec.el(coords_dof_idx);
        centroid[d] /= static_cast<double>(d_nodes_per_side);
    }

    // collect info across all processors
    const int total_num_elem_sides = SAMRAI_MPI::sumReduction(d_local_num_elem_sides);
    SAMRAI_MPI::sumReduction(&centroid[0], NDIM);

    // finish computing centroid
    for (int d = 0; d < NDIM; ++d)
    {
        centroid[d] /= static_cast<double>(total_num_elem_sides);
        d_current_location[d] = centroid[d];
    }
    return;
}

void
Source::updateSourceStrength(double loop_time, double dt, int source_number)
{
    // Q satisfies:
    //
    //    R Q + L dQ/dt = P_rsvr - P
    //
    // which is discretized as:
    //
    //    R (Q(n+1)+Q(n))/2 + L (Q(n+1)-Q(n))/dt = P_RSVR - P(n+1/2).
    double Q_src = d_q_src;
    double P_RSVR = 0.0;
    d_time = loop_time;
    if (loop_time < d_ramp_time && d_ramp_time > 0)
        P_RSVR = (loop_time / d_ramp_time) * d_P_reservoir;
    else
        P_RSVR = d_P_reservoir;
    
    // source is either a "pressure source" or a "flow source"
    if (d_is_pressure_source) 
    {
        // pressure source
        d_q_src = (2.0 * P_RSVR * dt - 2.0 * d_p_src * dt - d_R_src * Q_src * dt + 2.0 * d_L_src * Q_src) /
              (d_R_src * dt + 2.0 * d_L_src);
    }
    else 
    {
        // flow source
        double Q_RSVR = 0.0;
        if (loop_time < d_ramp_time && d_ramp_time > 0)
            Q_RSVR = (loop_time / d_ramp_time) * d_prescribed_flow;
        else
            Q_RSVR = d_prescribed_flow;
        
        d_q_src = Q_RSVR;
    }
   
    if (d_is_pressure_meter)
    {
        // zero out flow if all we care about is sampling the pressure
        d_q_src = 0.0;
    }
    
    pout << "Q_src[" << source_number << "] = " << setprecision(6) << setw(8) << d_q_src << " P_src[" << source_number
         << "] = " << setprecision(6) << setw(8) << d_p_src << " R_src[" << source_number << "] = " << setprecision(6)
         << setw(8) << d_R_src << " L_src[" << source_number << "] = " << setprecision(6) << setw(8) << d_L_src << "\n";

    // Write the current state to disk.
    writeDataFile();
    return;
}

Source::~Source()
{
}

void
Source::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_time", d_time);
    db->putDouble("d_q_src", d_q_src);
    db->putDouble("d_p_src", d_p_src);
    db->putDoubleArray("d_current_location", &d_current_location[0], NDIM);
    return;
} // putToDatabase

void
Source::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_time = db->getDouble("d_time");
    d_q_src = db->getDouble("d_q_src");
    d_p_src = db->getDouble("d_p_src");
    db->getDoubleArray("d_current_location", &d_current_location[0], NDIM);

    return;
} // getFromRestart

void
Source::writeDataFile()
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        d_data_stream.open(d_data_directory_name + "/" + d_data_filename, ios::app);
        d_data_stream.setf(ios_base::scientific);
        d_data_stream.precision(5);
        d_data_stream << d_time << " ";
        d_data_stream << d_q_src / flconv << " ";
        d_data_stream << d_p_src / prconv << "\n";
        d_data_stream.close();
    }
    return;
}
