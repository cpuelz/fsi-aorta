// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// APPLICATION INCLUDES
#include <four_chambered_heart/BoundaryConditions.h>
#include <four_chambered_heart/MechanicsModel.h>
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/ModelInitialization.h>
#include <four_chambered_heart/ModelParameters.h>
#include <four_chambered_heart/PartContext.h>

// BOOST INCLUDES
#include <boost/math/tools/roots.hpp>

// LIBMESH INCLUDES
#include <libmesh/boundary_info.h>
#include <libmesh/bounding_box.h>
#include <libmesh/centroid_partitioner.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/id_types.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/point.h>

// IBAMR includes
#include <ibtk/AppInitializer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <sstream>
#include <stddef.h>
#include <string>
#include <vector>
#include <set>

// initialize static map for active strain
ModelInitialization::ActiveStrainMap ModelInitialization::active_strain_model_map = { {"transverse",  Active_Strain_Model::TRANSVERSELY_ISOTROPIC},
                                                                                      {"orthotropic", Active_Strain_Model::ORTHOTROPIC},
                                                                                      {"transmural",  Active_Strain_Model::TRANSMURALLY_ORTHOTROPIC} };
// CLASS IMPLEMENTATION

void
ModelInitialization::transform_meshes(std::vector<MeshBase*>& mesh_vector,
                                      Pointer<AppInitializer> app_initializer,
                                      const std::vector<std::string>& part_names)
{
    pout << "\n";
    pout << "**************************************************************************** \n";
    pout << "transforming all the meshes...." << std::endl;
    pout << "**************************************************************************** \n";

    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // loop over the meshes corresponding to each part
    libMesh::Point lower(1e6, 1e6, 1e6);
    libMesh::Point upper(-1e6, -1e6, -1e6);
    for (unsigned int part = 0; part < part_names.size(); ++part)
    {
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        pout << "\n";
        pout << "mesh info for part " << part << "....\n";
        pout << "number of elems = " << mesh_vector[part]->n_elem() << "\n"
             << "number of nodes = " << mesh_vector[part]->n_nodes() << "\n"
             << "max elem ID = " << mesh_vector[part]->max_elem_id() << "\n"
             << "max node ID = " << mesh_vector[part]->max_node_id() << "\n";
        plog << "processor " << SAMRAI_MPI::getRank() << ": number of        local elems = "
             << distance(mesh_vector[part]->local_elements_begin(), mesh_vector[part]->local_elements_end()) << "\n";
        plog << "processor " << SAMRAI_MPI::getRank() << ": number of active local elems = "
             << distance(mesh_vector[part]->active_local_elements_begin(),
                         mesh_vector[part]->active_local_elements_end())
             << "\n";
        plog << "processor " << SAMRAI_MPI::getRank() << ": number of        local nodes = "
             << distance(mesh_vector[part]->local_nodes_begin(), mesh_vector[part]->local_nodes_end()) << "\n";

        // Rescale the mesh so that length is in cm.
        const string mesh_length_unit = part_db->getString("MESH_LENGTH_UNIT");
        pout << "input mesh length unit is: " << mesh_length_unit << "\n";
        if (strcasecmp(mesh_length_unit.c_str(), "m") == 0)
        {
            MeshTools::Modification::scale(*mesh_vector[part], 100.0);
            for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
            {
                if (MeshInfo::afterload_parms[ii].part == part)
                {
                    MeshInfo::afterload_parms[ii].radius *= 100.0;
                    MeshInfo::afterload_parms[ii].centroid *= 100.0;
                }
            }
        }
        else if (strcasecmp(mesh_length_unit.c_str(), "cm") == 0)
        {
            // do nothing
        }
        else if (strcasecmp(mesh_length_unit.c_str(), "mm") == 0)
        {
            MeshTools::Modification::scale(*mesh_vector[part], 0.1);
            for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
            {
                if (MeshInfo::afterload_parms[ii].part == part)
                {
                    MeshInfo::afterload_parms[ii].radius *= 0.1;
                    MeshInfo::afterload_parms[ii].centroid *= 0.1;
                }
            }
        }
        else if (strcasecmp(mesh_length_unit.c_str(), "1e-4 cm") == 0)
        {
            MeshTools::Modification::scale(*mesh_vector[part], 1e-4);
            for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
            {
                if (MeshInfo::afterload_parms[ii].part == part)
                {
                    MeshInfo::afterload_parms[ii].radius *= 1e-4;
                    MeshInfo::afterload_parms[ii].centroid *= 1e-4;
                }
            }
        }
        else
        {
            TBOX_ERROR("MESH_LENGTH_UNIT = " << mesh_length_unit << " is unrecognized" << endl);
        }

        // Center the meshes in the computational domain.
#if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
        const auto bbox = MeshTools::bounding_box(*mesh_vector[part]);
#else
        const auto bbox = MeshTools::create_bounding_box(*mesh_vector[part]);
#endif
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            lower(d) = min(bbox.min()(d), lower(d));
            upper(d) = max(bbox.max()(d), upper(d));
        }
    }

    pout << "\n";
    pout << "mesh bounding box lower before translation = (" << lower(0) << " , " << lower(1) << " , " << lower(2)
         << ") cm\n"
         << "mesh bounding box upper before translation = (" << upper(0) << " , " << upper(1) << " , " << upper(2)
         << ") cm\n";

    // translate all the meshes
    libMesh::Point lower_new(1e6, 1e6, 1e6);
    libMesh::Point upper_new(-1e6, -1e6, -1e6);
    for (unsigned int part = 0; part < mesh_vector.size(); ++part)
    {
        MeshTools::Modification::translate(*mesh_vector[part],
                                           -0.5 * (upper(0) + lower(0)),
                                           -0.5 * (upper(1) + lower(1)),
                                           -0.5 * (upper(2) + lower(2)));

#if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
        const auto bbox = MeshTools::bounding_box(*mesh_vector[part]);
#else
        const auto bbox = MeshTools::create_bounding_box(*mesh_vector[part]);
#endif
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            lower_new(d) = min(bbox.min()(d), lower_new(d));
            upper_new(d) = max(bbox.max()(d), upper_new(d));
        }
    }

    pout << "\n";
    pout << "mesh bounding box lower after translation = (" << lower_new(0) << " , " << lower_new(1) << " , "
         << lower_new(2) << ") cm\n"
         << "mesh bounding box upper after translation = (" << upper_new(0) << " , " << upper_new(1) << " , "
         << upper_new(2) << ") cm\n";

    // translate afterload model locations
    for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
    {
        MeshInfo::afterload_parms[ii].centroid(0) += -0.5 * (upper(0) + lower(0));
        MeshInfo::afterload_parms[ii].centroid(1) += -0.5 * (upper(1) + lower(1));
        MeshInfo::afterload_parms[ii].centroid(2) += -0.5 * (upper(2) + lower(2));
    }

    // Below is another translation to be sure that the extensions are flush with
    // the computational domain. Not sure how to do this in generality... this
    // code is very much geometry dependent.
    const double L = input_db->getDouble("L");
    libMesh::Point translation(-upper_new(0) + 9.0*L, 0.0, 0.0);
    for (unsigned int part = 0; part < mesh_vector.size(); ++part)
    {
        if (MeshInfo::afterload_parms.size() > 0)
        {
            MeshTools::Modification::translate(*mesh_vector[part], translation(0), translation(1), translation(2));
        }
    }
    for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
    {
        MeshInfo::afterload_parms[ii].centroid(0) += translation(0);
        MeshInfo::afterload_parms[ii].centroid(1) += translation(1);
        MeshInfo::afterload_parms[ii].centroid(2) += translation(2);
    }

    return;
}

// creates an equation systems object to store all the fiber information.
void
ModelInitialization::init_fiber_system(ExodusII_IO& mesh_reader,
                                       EquationSystems& system,
                                       Pointer<Database> input_db,
                                       const std::size_t sys_size)
{
    // initialize but do not populate equation system object to store fiber information
    const std::vector<std::string>& fiber_variable_names = mesh_reader.get_elem_var_names();
    if (fiber_variable_names.size() == 0)
        TBOX_ERROR(
            "ModelInitialization::init_fiber_system "
            ": mesh part doesn't have the fiber systems set up correctly, or they don't exist. This could be"
            " an issue with how the fibers were written to the exodus file.");
    if (fiber_variable_names.size() != sys_size)
        TBOX_ERROR(
            "ModelInitialization::init_fiber_system: "
            "mesh part has the incorrect system_data object size "
            << sys_size
            << " for the number of "
               "fiber systems "
            << fiber_variable_names.size() << ".");
    ExplicitSystem& fiber_info = system.add_system<ExplicitSystem>(input_db->getString("FIBER_SYSTEM_NAME"));
    pout << "fiber system variable names are: \n";
    for (unsigned int ii = 0; ii < fiber_variable_names.size(); ii++)
    {
        fiber_info.add_variable(fiber_variable_names[ii], CONSTANT, MONOMIAL);
        pout << fiber_variable_names[ii] << "\n";
    }
    return;
}

void
ModelInitialization::reinit_equation_systems(std::vector<EquationSystems*>& eq_system_vector,
                                             Pointer<AppInitializer> app_initializer,
                                             const std::vector<std::string>& part_names)
{
    for (unsigned int part = 0; part < part_names.size(); ++part)
    {
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        if (part_db->getBoolWithDefault("locally_refine_mesh", false))
        {
            pout << "\n";
            pout << "**************************************************************************** \n";
            eq_system_vector[part]->reinit();
            pout << "equation systems for " << part_names[part] << " have been reinitialized\n";
            pout << "**************************************************************************** \n";
        }
        if (part_db->getBoolWithDefault("uniformly_refine_mesh", false))
        {
            pout << "\n";
            pout << "**************************************************************************** \n";
            eq_system_vector[part]->reinit();
            pout << "equation systems for " << part_names[part] << " have been reinitialized\n";
            pout << "**************************************************************************** \n";
        }
    }
}

void
ModelInitialization::uniformly_refine_meshes(std::vector<MeshBase*>& mesh_vector,
                                             Pointer<AppInitializer> app_initializer,
                                             const std::vector<std::string>& part_names)
{
    for (unsigned int part = 0; part < part_names.size(); ++part)
    {
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        if (part_db->getBoolWithDefault("uniformly_refine_mesh", false))
        {
            pout << "\n";
            pout << "**************************************************************************** \n";
            pout << "uniformly refining mesh part " << part << std::endl;
            MeshRefinement(*mesh_vector[part]).uniformly_refine();
            pout << "mesh for " << part_names[part] << " has been uniformly refined"
                 << "\n";
            pout << "**************************************************************************** \n";
        }
    }
}

void
ModelInitialization::locally_refine_meshes(std::vector<MeshBase*>& mesh_vector,
                                             Pointer<AppInitializer> app_initializer,
                                             const std::vector<std::string>& part_names)
{
    for (unsigned int part = 0; part < part_names.size(); ++part)
    {
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        if (part_db->getBoolWithDefault("locally_refine_mesh", false))
        {
            std::set<libMesh::subdomain_id_type> refinement_ids(MeshInfo::local_refinement_subdomain_IDs[part].begin(),
                                                                MeshInfo::local_refinement_subdomain_IDs[part].end());
            pout << "\n";
            pout << "**************************************************************************** \n";
            pout << "locally refining mesh part " << part << std::endl;
            libMesh::MeshRefinement mesh_refinement(*mesh_vector[part]);
            // libMesh complains if flags are set inconsistently between
            // processors, so just tag all cells we want to refine on all
            // processors:
            libMesh::MeshBase::const_element_iterator el = mesh_vector[part]->active_elements_begin();
            const libMesh::MeshBase::const_element_iterator end_el = mesh_vector[part]->active_elements_end();
            for (; el != end_el; ++el)
            {
                libMesh::Elem * elem = *el;
                const libMesh::subdomain_id_type block_id = elem->subdomain_id();
                if (refinement_ids.count(block_id) != 0)
                {
                    elem->set_refinement_flag(libMesh::Elem::REFINE);
                }
                else
                {
                    elem->set_refinement_flag(libMesh::Elem::DO_NOTHING);
                }
            }
            // refine mesh
            mesh_refinement.refine_elements();
            pout << "mesh for " << part_names[part] << " has been locally refined"
                    << "\n";
            pout << "**************************************************************************** \n";
        }
    }
}

// for exodus meshes
void
ModelInitialization::load_fiber_system(ExodusII_IO& mesh_reader, EquationSystems& system, Pointer<Database> input_db)
{
    // populate equation system object to store fiber information
    const std::vector<std::string>& fiber_variable_names = mesh_reader.get_elem_var_names();
    ExplicitSystem& fiber_info = system.get_system<ExplicitSystem>(input_db->getString("FIBER_SYSTEM_NAME"));
    for (unsigned int ii = 0; ii < fiber_variable_names.size(); ii++)
    {
        mesh_reader.copy_elemental_solution(fiber_info, fiber_variable_names[ii], fiber_variable_names[ii]);
    }

    return;
}

// this function populates the MeshInfo object with blocks, sidesets, and nodesets, i.e.
// their IDs and corresponding names.
void
ModelInitialization::load_mesh_info(std::vector<MeshBase*>& mesh_vector,
                                    Pointer<AppInitializer> app_initializer,
                                    const std::vector<std::string>& part_names)
{
    MeshInfo::has_fibers.resize(mesh_vector.size());
    MeshInfo::turn_on_meters.resize(mesh_vector.size());
    MeshInfo::turn_on_sources.resize(mesh_vector.size());
    MeshInfo::subdomain_names.resize(mesh_vector.size());
    MeshInfo::circ_model_body_force_IDs.resize(mesh_vector.size());
    MeshInfo::subdomain_IDs.resize(mesh_vector.size());
    Pointer<Database> afterload_db = app_initializer->getComponentDatabase("AfterloadModel");
    Pointer<Database> model_parm_db = app_initializer->getComponentDatabase("ModelParameters");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    const double DX = input_db->getDouble("DX");

    // loop over the meshes corresponding to each part
    for (int part = 0; part < static_cast<int>(mesh_vector.size()); ++part)
    {
        // figure out if this mesh part has fibers, meters, sources.
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        MeshInfo::has_fibers[part] = part_db->getBoolWithDefault("has_fibers", false);
        MeshInfo::turn_on_meters[part] = part_db->getBoolWithDefault("turn_on_meters", false);
        MeshInfo::turn_on_sources[part] = part_db->getBoolWithDefault("turn_on_sources", false);

        // get block ids and names from mesh file
        std::map<short unsigned int, std::string> temp_block_map = mesh_vector[part]->get_subdomain_name_map();
        std::set<subdomain_id_type> set_of_subdomain_ids;
        mesh_vector[part]->subdomain_ids(set_of_subdomain_ids);
        std::map<std::string, int> temp_map;

        // first loop over the block map and see if all subdomains actually have string names
        for (std::set<subdomain_id_type>::iterator it = set_of_subdomain_ids.begin(); it != set_of_subdomain_ids.end();
             ++it)
        {
            const std::string subdomain_name = temp_block_map[*it];
            if (subdomain_name.empty())
            {
                std::stringstream it_string;
                it_string << *it;
                temp_block_map[*it] = "default_subdomain_name_" + it_string.str();
                TBOX_WARNING(
                    "ModelInitialization::load_mesh_info: "
                    "Subdomain "
                    << *it << " in part " << part
                    << " is not named. "
                       "You need to add a database for this subdomain called: "
                    << temp_block_map[*it] << " in the input file.");
            }
        }

        // swap keys and values, and store in MeshInfo object
        std::map<short unsigned int, std::string>::iterator bb;
        for (bb = temp_block_map.begin(); bb != temp_block_map.end(); ++bb)
        {
            temp_map.insert(std::pair<std::string, int>(bb->second, (int)bb->first));
            MeshInfo::subdomain_names[part].push_back(bb->second);
            MeshInfo::subdomain_IDs[part].push_back(bb->first);
        }
        MeshInfo::block_map.push_back(temp_map);
        temp_map.clear();

        // get sideset information *******
        BoundaryConditions::boundary_info.push_back(&(mesh_vector[part]->get_boundary_info()));

        std::map<boundary_id_type, std::string> temp_sideset_map =
            mesh_vector[part]->get_boundary_info().get_sideset_name_map();

        // first loop over the sideset map and see if all sidesets actually have string names
        const BoundaryInfo& boundary_info = mesh_vector[part]->get_boundary_info();
        for (std::set<boundary_id_type>::iterator it = boundary_info.get_side_boundary_ids().begin(); it != boundary_info.get_side_boundary_ids().end();
             ++it)
        {
            const std::string sideset_name = boundary_info.get_sideset_name(*it);
            if (sideset_name.empty())
            {
                std::stringstream it_string;
                it_string << *it;
                temp_sideset_map[*it] = "default_sideset_name_" + it_string.str();
                TBOX_WARNING(
                    "ModelInitialization::load_mesh_info: "
                    "Sideset "
                    << *it << " in part " << part
                    << " is not named. Naming it " << temp_sideset_map[*it]);
            }
        }

        // swap keys and values, and store in MeshInfo object
        std::map<boundary_id_type, std::string>::iterator ss;
        for (ss = temp_sideset_map.begin(); ss != temp_sideset_map.end(); ++ss)
        {
            temp_map.insert(std::pair<std::string, int>(ss->second, (int)ss->first));
        }
        MeshInfo::sideset_map.push_back(temp_map);
        temp_map.clear();

        // get nodeset information *******
        std::map<boundary_id_type, std::string> temp_nodeset_map =
            mesh_vector[part]->get_boundary_info().get_nodeset_name_map();

        // check to see if specified nodeset IDs for the meters actually exist in the mesh
        if (MeshInfo::turn_on_meters[part])
        {
            SAMRAI::tbox::Array<int> meter_IDs =  part_db->getIntegerArray("nodeset_IDs_for_meters");
            for(int mm = 0; mm < meter_IDs.getSize(); ++mm)
            {
                if(boundary_info.get_node_boundary_ids().find(meter_IDs[mm]) == boundary_info.get_node_boundary_ids().end())
                {
                    TBOX_ERROR("ModelInitialization::load_mesh_info: nodeset ID " << meter_IDs[mm]
                            << " corresponding to a meter was not found in list of nodesets IDs from the mesh."
                            << " check to make sure you have the correct nodeset ID specified in the input file.");
                }
            }
        }

        // first loop over the nodeset map and see if all nodesets actually have string names
        for (std::set<boundary_id_type>::iterator it = boundary_info.get_node_boundary_ids().begin(); it != boundary_info.get_node_boundary_ids().end();
             ++it)
        {
            const std::string nodeset_name = boundary_info.get_nodeset_name(*it);
            if (nodeset_name.empty())
            {
                std::stringstream it_string;
                it_string << *it;
                temp_nodeset_map[*it] = "default_nodeset_name_" + it_string.str();
                TBOX_WARNING(
                    "ModelInitialization::load_mesh_info: "
                    "Nodeset "
                    << *it << " in part " << part
                    << " is not named. Naming it " << temp_nodeset_map[*it]);
            }
        }

        // swap keys and values, and store in MeshInfo object
        std::map<boundary_id_type, std::string>::iterator nn;
        for (nn = temp_nodeset_map.begin(); nn != temp_nodeset_map.end(); ++nn)
        {
            temp_map.insert(std::pair<std::string, int>(nn->second, (int)nn->first));
        }
        MeshInfo::nodeset_map.push_back(temp_map);

        pout << "\n";
        pout << "**************************************************************************** \n";
        pout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART = " << part << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
        pout << "**************************************************************************** \n";

        pout << "\n";
        pout << "Nodeset info for part " << part << "...." << std::endl;
        std::map<std::string, int>::iterator n;
        for (n = MeshInfo::nodeset_map[part].begin(); n != MeshInfo::nodeset_map[part].end(); ++n)
        {
            pout << n->first << " " << n->second << std::endl;
        }

        pout << "\n";
        pout << "Sideset info for part " << part << "...." << std::endl;
        std::map<std::string, int>::iterator s;
        for (s = MeshInfo::sideset_map[part].begin(); s != MeshInfo::sideset_map[part].end(); ++s)
        {
            pout << s->first << " " << s->second << std::endl;
        }

        pout << "\n";
        pout << "Block info for part " << part << "...." << std::endl;
        std::map<std::string, int>::iterator b;
        for (b = MeshInfo::block_map[part].begin(); b != MeshInfo::block_map[part].end(); ++b)
        {
            pout << b->first << " " << b->second << std::endl;
        }

        MeshInfo::tether_body_force_names.push_back(part_db->getStringArray("tether_body_force_names"));
        MeshInfo::tether_surface_force_names.push_back(part_db->getStringArray("tether_surface_force_names"));
        MeshInfo::surface_pressure_names.push_back(part_db->getStringArray("surface_pressure_names"));
        if (part_db->getBoolWithDefault("use_afterload_models", false))
        {
            MeshInfo::afterload_model_nodeset_names.push_back(part_db->getStringArray("afterload_model_nodeset_names"));
        }
        else
        {
            SAMRAI::tbox::Array<std::string> empty_string_vector;
            MeshInfo::afterload_model_nodeset_names.push_back(empty_string_vector);
        }
        if (part_db->getBoolWithDefault("locally_refine_mesh", false))
        {
            MeshInfo::local_refinement_subdomain_names.push_back(part_db->getStringArray("local_refinement_subdomain_names"));
        }
        else
        {
            SAMRAI::tbox::Array<std::string> empty_string_vector;
            MeshInfo::local_refinement_subdomain_names.push_back(empty_string_vector);
        }

        // for body tether forces
        std::vector<int> temp;
        for (int ii = 0; ii < MeshInfo::tether_body_force_names[part].getSize(); ++ii)
        {
            temp.push_back(MeshInfo::block_map[part][MeshInfo::tether_body_force_names[part][ii]]);
        }
        MeshInfo::tether_body_force_IDs.push_back(temp);
        temp.clear();

        // for surface tether forces
        for (int ii = 0; ii < MeshInfo::tether_surface_force_names[part].getSize(); ++ii)
        {
            temp.push_back(MeshInfo::sideset_map[part][MeshInfo::tether_surface_force_names[part][ii]]);
        }
        MeshInfo::tether_surface_force_IDs.push_back(temp);
        temp.clear();

        // for surface pressure
        for (int ii = 0; ii < MeshInfo::surface_pressure_names[part].getSize(); ++ii)
        {
            temp.push_back(MeshInfo::sideset_map[part][MeshInfo::surface_pressure_names[part][ii]]);
        }
        MeshInfo::surface_pressure_IDs.push_back(temp);

        // for afterload models
        if (part_db->getBoolWithDefault("use_afterload_models", false))
        {
            MeshInfo::circ_model_surface_tethering = afterload_db->getBool("surface_tethering");
            MeshInfo::circ_model_body_tethering = afterload_db->getBool("body_tethering");
            MeshInfo::Cartesian_L[0] = afterload_db->getDouble("Cartesian_L_x");
            MeshInfo::Cartesian_L[1] = afterload_db->getDouble("Cartesian_L_y");
            MeshInfo::Cartesian_L[2] = afterload_db->getDouble("Cartesian_L_z");
            for (int ii = 0; ii < MeshInfo::afterload_model_nodeset_names[part].getSize(); ++ii)
            {
                Afterload_Parms afterload_parms;
                afterload_parms.nodeset_name = MeshInfo::afterload_model_nodeset_names[part][ii];
                Pointer<Database> src_db = afterload_db->getDatabase(afterload_parms.nodeset_name);

                // get the ID for the nodeset at the top of the PA.
                afterload_parms.nodeset_ID = MeshInfo::nodeset_map[part][afterload_parms.nodeset_name];

                // compute the centroid for the afterload model
                int count = 0;
                libMesh::Point sum(0.0, 0.0, 0.0);
                for (unsigned int mm = 0; mm < mesh_vector[part]->n_nodes(); ++mm)
                {
                    const Node* node_ptr = mesh_vector[part]->node_ptr(mm);
                    std::vector<short int> bdry_ids;
                    BoundaryConditions::boundary_info[part]->boundary_ids(node_ptr, bdry_ids);
                    if (find(bdry_ids.begin(), bdry_ids.end(), afterload_parms.nodeset_ID) != bdry_ids.end())
                    {
                        count += 1;
                        sum += *node_ptr;
                    }
                }
                sum /= static_cast<double>(count);
                afterload_parms.centroid = sum;

                // compute the radius for the afterload model
                double radius = std::numeric_limits<double>::max();
                for (unsigned int mm = 0; mm < mesh_vector[part]->n_nodes(); ++mm)
                {
                    const Node* node_ptr = mesh_vector[part]->node_ptr(mm);
                    std::vector<short int> bdry_ids;
                    BoundaryConditions::boundary_info[part]->boundary_ids(node_ptr, bdry_ids);
                    if (find(bdry_ids.begin(), bdry_ids.end(), afterload_parms.nodeset_ID) != bdry_ids.end())
                    {
                        radius = min(radius, (*node_ptr - sum).norm());
                    }
                }
		const double buffer = src_db->getDouble("radius_buffer");
                afterload_parms.radius = radius - buffer * DX;
                afterload_parms.axis = src_db->getInteger("axis");
                afterload_parms.side = src_db->getInteger("side");
                afterload_parms.part = part;
                const std::string subdomain_name = src_db->getString("subdomain_name");
                MeshInfo::circ_model_body_force_IDs[part].push_back(MeshInfo::block_map[part][subdomain_name]);
                const int subdomain_ID = MeshInfo::block_map[part][subdomain_name];
                MeshInfo::subdomain_to_afterload_parms.insert(std::pair<int, Afterload_Parms>(subdomain_ID, afterload_parms));
                MeshInfo::afterload_parms.push_back(afterload_parms);
            }
        }
        // for local mesh refinement
        if (part_db->getBoolWithDefault("locally_refine_mesh", false))
        {
            std::vector<int> temp;
            for (int ii = 0; ii < MeshInfo::local_refinement_subdomain_names[part].getSize(); ++ii)
            {
                temp.push_back(MeshInfo::block_map[part][MeshInfo::local_refinement_subdomain_names[part][ii]]);
            }
            MeshInfo::local_refinement_subdomain_IDs.push_back(temp);
            temp.clear();
        }

        // print out stuff.
        pout << "\n";
        pout << "body tether force IDs for part " << part << "...." << std::endl;
        if (MeshInfo::tether_body_force_IDs[part][0] == 0)
        {
            pout << "NONE\n";
        }
        else
        {
            for (unsigned int ii = 0; ii < MeshInfo::tether_body_force_IDs[part].size(); ++ii)
            {
                pout << MeshInfo::tether_body_force_IDs[part][ii] << " ( "
                     << MeshInfo::tether_body_force_names[part][ii] << " ) "
                     << "\n";
            }
        }

        pout << "\n";
        pout << "surface tether force IDs for part " << part << "...." << std::endl;
        if (MeshInfo::tether_surface_force_IDs[part][0] == 0)
        {
            pout << "NONE\n";
        }
        else
        {
            for (unsigned int ii = 0; ii < MeshInfo::tether_surface_force_IDs[part].size(); ++ii)
            {
                pout << MeshInfo::tether_surface_force_IDs[part][ii] << " ( "
                     << MeshInfo::tether_surface_force_names[part][ii] << " ) "
                     << "\n";
            }
        }

        if (part_db->getBoolWithDefault("locally_refine_mesh", false))
        {
            pout << "\n";
            pout << "local mesh refinement IDs for part " << part << "...." << std::endl;
            for (unsigned int ii = 0; ii < MeshInfo::local_refinement_subdomain_IDs[part].size(); ++ii)
            {
                pout << MeshInfo::local_refinement_subdomain_IDs[part][ii] << " ( "
                        << MeshInfo::local_refinement_subdomain_names[part][ii] << " ) "
                        << "\n";
            }
        }


	MeshInfo::pericardial_tethering_names = model_parm_db->getStringArray("pericardial_tethering_names");
	// for pericardial tethering
	for (int ii = 0; ii < MeshInfo::pericardial_tethering_names.getSize(); ++ii)
	  {
	    MeshInfo::pericardial_tethering_IDs.push_back(
							  MeshInfo::sideset_map[part][MeshInfo::pericardial_tethering_names[ii]]);
	  }
	
	pout << "\n";
	pout << "pericardial tethering IDs for the part " << part << "...." << std::endl;
	if (MeshInfo::pericardial_tethering_IDs[0] == 0)
	  {
	    pout << "NONE\n";
	  }
	else
	  {
	    for (unsigned int ii = 0; ii < MeshInfo::pericardial_tethering_IDs.size(); ++ii)
	      {
		pout << MeshInfo::pericardial_tethering_IDs[ii] << " ( "
		     << MeshInfo::pericardial_tethering_names[ii] << " ) "
		     << "\n";
	      }
	  }
        
    } // end loop over parts

    pout << "\n";
    pout << "**************************************************************************** \n";
    pout << "afterload models are...." << std::endl;
    if (MeshInfo::afterload_parms.size() == 0)
    {
        pout << "NONE\n";
    }
    else
    {
        for (unsigned int ii = 0; ii < MeshInfo::afterload_parms.size(); ++ii)
        {
            pout << "name: " << MeshInfo::afterload_parms[ii].nodeset_name << "\n";
            pout << "part: " << MeshInfo::afterload_parms[ii].part << "\n";
            pout << "axis: " << MeshInfo::afterload_parms[ii].axis << "\n";
            pout << "side: " << MeshInfo::afterload_parms[ii].side << "\n";
        }
    }
    pout << "**************************************************************************** \n";
    pout << "\n";

    return;
}

// function to load model parameters
void
ModelInitialization::load_model_parameters(Pointer<AppInitializer> app_initializer,
                                           const std::vector<std::string>& part_names,
                                           std::vector<std::unique_ptr<PartContext> >& context_vector)
{
    // loop over parts
    for (unsigned int part = 0; part < part_names.size(); ++part)
    {
        Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
        // loop over subdomains
        for (int ss = 0; ss < static_cast<int>(MeshInfo::subdomain_names[part].size()); ++ss)
        {
            const std::string subdomain_name = MeshInfo::subdomain_names[part][ss];
            const libMesh::subdomain_id_type subdomain_ID = MeshInfo::block_map[part][subdomain_name];
            Pointer<Database> subdomain_db = part_db->getDatabase(subdomain_name);

            std::unique_ptr<General_ModelParms> general_ptr(new General_ModelParms);
            std::unique_ptr<Guccione_ModelParms> guccione_ptr(new Guccione_ModelParms);
            std::unique_ptr<HGO_ModelParms> hgo_ptr(new HGO_ModelParms);
            std::unique_ptr<HO_ModelParms> ho_ptr(new HO_ModelParms);
            std::unique_ptr<Neohookean_ModelParms> neohookean_ptr(new Neohookean_ModelParms);
            std::unique_ptr<Exp_Neohookean_ModelParms> exp_neohookean_ptr(new Exp_Neohookean_ModelParms);
            std::unique_ptr<Biaxial_Fung_Type_ModelParms> biaxial_fung_type_ptr(new Biaxial_Fung_Type_ModelParms);
            std::unique_ptr<Augustin_ModelParms> augustin_ptr(new Augustin_ModelParms);
            std::unique_ptr<NonlinearSpring_ModelParms> nonlinear_spring_ptr(new NonlinearSpring_ModelParms);

            // setup stuff for logging invariants during the simulation
            context_vector[part]->min_J[subdomain_ID] = 0.5 * std::numeric_limits<double>::max();
            context_vector[part]->max_J[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();
            context_vector[part]->min_I1[subdomain_ID] = 0.5 * std::numeric_limits<double>::max();
            context_vector[part]->max_I1[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();
            context_vector[part]->min_I4_one[subdomain_ID] = 0.5 * std::numeric_limits<double>::max();
            context_vector[part]->max_I4_one[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();
            context_vector[part]->min_I4_two[subdomain_ID] = 0.5 * std::numeric_limits<double>::max();
            context_vector[part]->max_I4_two[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();
            context_vector[part]->min_I4_three[subdomain_ID] = 0.5 * std::numeric_limits<double>::max();
            context_vector[part]->max_I4_three[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();

            // setup stuff for logging active tension
            context_vector[part]->active_tension[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();

            // setting up the dilational stress type.
            general_ptr->beta_s = subdomain_db->getDoubleWithDefault("beta_s", 1.0);
            if (subdomain_db->getString("volumetric_energy_type")
                    .compare(MechanicsModel::LOG_VOLUMETRIC_ENERGY_NAME))
            {
                context_vector[part]->volumetric_energy_ID[subdomain_ID] = MechanicsModel::LOG_VOLUMETRIC_ENERGY_ID;
            }
            else if (subdomain_db->getString("volumetric_energy_type")
                    .compare(MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_NAME))
            {
                context_vector[part]->volumetric_energy_ID[subdomain_ID] = MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_ID;
            }
            else
            {
                TBOX_ERROR("ModelInitialization::load_model_parameters: volumetric energy type not recognized.");
            }

            general_ptr->enable_active_stress = subdomain_db->getBoolWithDefault("enable_active_stress", false);
            if (general_ptr->enable_active_stress)
            {
                general_ptr->t_ramp = subdomain_db->getDouble("t_ramp");
                general_ptr->nu1 = subdomain_db->getDouble("nu1");
                general_ptr->nu2 = subdomain_db->getDouble("nu2");
                general_ptr->nu3 = subdomain_db->getDouble("nu3");
                general_ptr->scaling_with_I4 = subdomain_db->getBool("scaling_with_I4");
                general_ptr->Tension = subdomain_db->getDouble("Tension");
                general_ptr->first_active_tension_fiber_index = subdomain_db->getInteger("first_active_tension_fiber_index");
                general_ptr->second_active_tension_fiber_index = subdomain_db->getInteger("second_active_tension_fiber_index");
                general_ptr->third_active_tension_fiber_index = subdomain_db->getInteger("third_active_tension_fiber_index");
                if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::original_active_tension_function;
                    general_ptr->original_active_tension_parms.t_active = active_tension_db->getDouble("t_active");
                    general_ptr->original_active_tension_parms.t_relax = active_tension_db->getDouble("t_relax");
                    general_ptr->original_active_tension_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->original_active_tension_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->original_active_tension_parms.t_delay = active_tension_db->getDouble("t_delay");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::augustin_active_tension_function;
                    general_ptr->augustin_active_tension_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->augustin_active_tension_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->augustin_active_tension_parms.t_delay = active_tension_db->getDouble("t_delay");
                    general_ptr->augustin_active_tension_parms.t_t = active_tension_db->getDouble("t_t");
                    general_ptr->augustin_active_tension_parms.tau_c = active_tension_db->getDouble("tau_c");
                    general_ptr->augustin_active_tension_parms.tau_r = active_tension_db->getDouble("tau_r");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::DATA_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::DATA_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::data_activation_function;
                    general_ptr->data_activation_function_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->data_activation_function_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->data_activation_function_parms.t_delay = active_tension_db->getDouble("t_delay");
                    general_ptr->data_activation_function_parms.t_active = active_tension_db->getDouble("t_active");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::PRESSURE_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::PRESSURE_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::pressure_activation_function;
                    general_ptr->pressure_activation_function_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->pressure_activation_function_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->pressure_activation_function_parms.t_delay = active_tension_db->getDouble("t_delay");
                    general_ptr->pressure_activation_function_parms.t_peak = active_tension_db->getDouble("t_peak");
                    general_ptr->pressure_activation_function_parms.t_plateau = active_tension_db->getDouble("t_plateau");
                    general_ptr->pressure_activation_function_parms.t_drop = active_tension_db->getDouble("t_drop");
                }
                 else
                {
                    TBOX_ERROR("ModelInitialization::load_model_parameters: active tension function type not recognized.");
                }
            }

            // Check if we are using the active strain from the input file
            general_ptr->enable_active_strain = subdomain_db->getBoolWithDefault("enable_active_strain", false);
            if (general_ptr->enable_active_strain)
            {
                general_ptr->Tension = subdomain_db->getDouble("Tension");
                general_ptr->first_active_tension_fiber_index = subdomain_db->getInteger("first_active_tension_fiber_index");
                general_ptr->second_active_tension_fiber_index = subdomain_db->getInteger("second_active_tension_fiber_index");
                general_ptr->third_active_tension_fiber_index = subdomain_db->getInteger("third_active_tension_fiber_index");
                // This is used only for the orthotropic and transmural models
                general_ptr->active_strain_kappa = subdomain_db->getDouble("active_strain_kappa");
                // how much f0 is going to shorten
                general_ptr->max_fiber_shortening = subdomain_db->getDouble("max_fiber_shortening");
                // find out which active strain model we are using
                std::string active_strain_model = subdomain_db->getString("active_strain_model");
                auto it_active_strain =  active_strain_model_map.find(active_strain_model);
                // set it if we find it
                if( it_active_strain != active_strain_model_map.end() )
                {
                    general_ptr->active_strain_model = it_active_strain->second;
                }
                // throw an error otherwise
                else
                {
                    pout << "ERROR: Active Strain model enabled but not correctly specified!" << std::endl;
                    pout << "Choose between: " << std::endl;
                    for(auto && m : active_strain_model_map ) std::cout << m.first << std::endl;
                    pout << std::endl;
                    throw std::runtime_error("Active Strain model enabled but not correctly specified!");
                }
                // register activation functions
                if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::original_active_tension_function;
                    general_ptr->original_active_tension_parms.t_active = active_tension_db->getDouble("t_active");
                    general_ptr->original_active_tension_parms.t_relax = active_tension_db->getDouble("t_relax");
                    general_ptr->original_active_tension_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->original_active_tension_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->original_active_tension_parms.t_delay = active_tension_db->getDouble("t_delay");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::augustin_active_tension_function;
                    general_ptr->augustin_active_tension_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->augustin_active_tension_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->augustin_active_tension_parms.t_delay = active_tension_db->getDouble("t_delay");
                    general_ptr->augustin_active_tension_parms.t_t = active_tension_db->getDouble("t_t");
                    general_ptr->augustin_active_tension_parms.tau_c = active_tension_db->getDouble("tau_c");
                    general_ptr->augustin_active_tension_parms.tau_r = active_tension_db->getDouble("tau_r");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                        .compare(MechanicsModel::DATA_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                    Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::DATA_ACTIVE_TENSION_FUNCTION_NAME);
                    general_ptr->active_tension_function = &MechanicsModel::data_activation_function;
                    general_ptr->data_activation_function_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                    general_ptr->data_activation_function_parms.t_period = active_tension_db->getDouble("t_period");
                    general_ptr->data_activation_function_parms.t_delay = active_tension_db->getDouble("t_delay");
                    general_ptr->data_activation_function_parms.t_active = active_tension_db->getDouble("t_active");
                }
                else if (subdomain_db->getString("active_tension_function_type")
                    .compare(MechanicsModel::PRESSURE_ACTIVE_TENSION_FUNCTION_NAME)== 0)
                {
                Pointer<Database> active_tension_db = subdomain_db->getDatabase(MechanicsModel::PRESSURE_ACTIVE_TENSION_FUNCTION_NAME);
                general_ptr->active_tension_function = &MechanicsModel::pressure_activation_function;
                general_ptr->pressure_activation_function_parms.t_ramp = active_tension_db->getDouble("t_ramp");
                general_ptr->pressure_activation_function_parms.t_period = active_tension_db->getDouble("t_period");
                general_ptr->pressure_activation_function_parms.t_delay = active_tension_db->getDouble("t_delay");
                general_ptr->pressure_activation_function_parms.t_peak = active_tension_db->getDouble("t_peak");
                general_ptr->pressure_activation_function_parms.t_plateau = active_tension_db->getDouble("t_plateau");
                general_ptr->pressure_activation_function_parms.t_drop = active_tension_db->getDouble("t_drop");
               }
             }

            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::HO_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::HO_NAME);
                general_ptr->t_ramp = subdomain_db->getDouble("t_ramp");
                ho_ptr->mu_e = model_db->getDouble("mu_e");
                ho_ptr->a = model_db->getDouble("a");
                ho_ptr->b = model_db->getDouble("b");
                ho_ptr->af = model_db->getDouble("af");
                ho_ptr->bf = model_db->getDouble("bf");
                ho_ptr->as = model_db->getDouble("as");
                ho_ptr->bs = model_db->getDouble("bs");
                ho_ptr->afs = model_db->getDouble("afs");
                ho_ptr->bfs = model_db->getDouble("bfs");
                ho_ptr->kf = model_db->getDouble("kf");
                ho_ptr->ks = model_db->getDouble("ks");
                ho_ptr->enable_prestrain = model_db->getBool("enable_prestrain");
                ho_ptr->radial_prestrain_coeff = model_db->getDouble("radial_prestrain_coeff");
                ho_ptr->circ_prestrain_coeff = model_db->getDouble("circ_prestrain_coeff");
                ho_ptr->radial_fiber_index = model_db->getInteger("radial_fiber_index");
                ho_ptr->circ_fiber_index = model_db->getInteger("circ_fiber_index");
                ho_ptr->single_fiber_family = model_db->getBool("single_fiber_family");
                ho_ptr->use_scaled_invariants = model_db->getBool("use_scaled_invariants");
                ho_ptr->turn_off_fibers_in_compression = model_db->getBool("turn_off_fibers_in_compression");
                ho_ptr->DEV_projection = model_db->getBool("DEV_projection");
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::HO_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::HO_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::HGO_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::HGO_NAME);
                hgo_ptr->c = model_db->getDouble("c");
                hgo_ptr->c0 = model_db->getDouble("c0");
                hgo_ptr->c2 = model_db->getDouble("c2");
                hgo_ptr->radial_prestrain_coeff = model_db->getDouble("radial_prestrain_coeff");
                hgo_ptr->circ_prestrain_coeff = model_db->getDouble("circ_prestrain_coeff");
                hgo_ptr->radial_fiber_index = model_db->getInteger("radial_fiber_index");
                hgo_ptr->circ_fiber_index = model_db->getInteger("circ_fiber_index");
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::HGO_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::HGO_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::GUCCIONE_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::GUCCIONE_NAME);
                guccione_ptr->C = model_db->getDouble("C");
                guccione_ptr->bf = model_db->getDouble("bf");
                guccione_ptr->bt = model_db->getDouble("bt");
                guccione_ptr->bfs = model_db->getDouble("bfs");
                guccione_ptr->circ_prestrain_coeff = model_db->getDouble("circ_prestrain_coeff");
                guccione_ptr->circ_fiber_index = model_db->getInteger("circ_fiber_index");
                guccione_ptr->DEV_projection = model_db->getBoolWithDefault("DEV_projection", true);
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::GUCCIONE_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::GUCCIONE_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::NEOHOOKEAN_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::NEOHOOKEAN_NAME);
                neohookean_ptr->mu_s = model_db->getDouble("mu_s");
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::NEOHOOKEAN_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::NEOHOOKEAN_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::EXP_NEOHOOKEAN_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::EXP_NEOHOOKEAN_NAME);
                exp_neohookean_ptr->a = model_db->getDouble("a");
                exp_neohookean_ptr->b = model_db->getDouble("b");
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::EXP_NEOHOOKEAN_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::EXP_NEOHOOKEAN_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::BIAXIAL_FUNG_TYPE_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::BIAXIAL_FUNG_TYPE_NAME);
                biaxial_fung_type_ptr->C = model_db->getDouble("C");
                biaxial_fung_type_ptr->A1 = model_db->getDouble("A1");
                biaxial_fung_type_ptr->A2 = model_db->getDouble("A2");
                biaxial_fung_type_ptr->A3 = model_db->getDouble("A3");
                biaxial_fung_type_ptr->A4 = model_db->getDouble("A4");
                biaxial_fung_type_ptr->A5 = model_db->getDouble("A5");
                biaxial_fung_type_ptr->A6 = model_db->getDouble("A6");
                biaxial_fung_type_ptr->radial_fiber_index = model_db->getInteger("radial_fiber_index");
                biaxial_fung_type_ptr->circ_fiber_index = model_db->getInteger("circ_fiber_index");
                biaxial_fung_type_ptr->DEV_projection = model_db->getBoolWithDefault("DEV_projection", true);
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::BIAXIAL_FUNG_TYPE_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::BIAXIAL_FUNG_TYPE_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::AUGUSTIN_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::AUGUSTIN_NAME);
                augustin_ptr->kappa = model_db->getDouble("kappa");
                augustin_ptr->a = model_db->getDouble("a");
                augustin_ptr->b = model_db->getDouble("b");
                augustin_ptr->af = model_db->getDouble("af");
                augustin_ptr->bf = model_db->getDouble("bf");
                augustin_ptr->circ_fiber_index = model_db->getInteger("circ_fiber_index");
                augustin_ptr->DEV_projection = model_db->getBoolWithDefault("DEV_projection", true);
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::AUGUSTIN_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::AUGUSTIN_NAME;
            }
            if (subdomain_db->getString("constitutive_model").compare(MechanicsModel::NONLINEAR_SPRING_NAME) == 0)
            {
                Pointer<Database> model_db = subdomain_db->getDatabase(MechanicsModel::NONLINEAR_SPRING_NAME);
                nonlinear_spring_ptr->af = model_db->getDouble("af");
                nonlinear_spring_ptr->bf = model_db->getDouble("bf");
                nonlinear_spring_ptr->fiber_index = model_db->getInteger("fiber_index");
                nonlinear_spring_ptr->DEV_projection = model_db->getBoolWithDefault("DEV_projection", true);
                nonlinear_spring_ptr->use_scaled_invariants = model_db->getBool("use_scaled_invariants");
                if (nonlinear_spring_ptr->use_scaled_invariants) TBOX_ERROR("ModelInitialization::load_model_parameters: nonlinear spring model is NOT implemented with scaled invariants.");
                nonlinear_spring_ptr->prestrain_coeff = model_db->getDouble("prestrain_coeff");
                nonlinear_spring_ptr->enable_prestrain = model_db->getBool("enable_prestrain");
                context_vector[part]->constitutive_model_ID[subdomain_ID] = MechanicsModel::NONLINEAR_SPRING_ID;
                context_vector[part]->constitutive_model_name[subdomain_ID] = MechanicsModel::NONLINEAR_SPRING_NAME;
            }

            context_vector[part]->general_parms[subdomain_ID] = std::move(general_ptr);
            context_vector[part]->guccione_parms[subdomain_ID] = std::move(guccione_ptr);
            context_vector[part]->ho_parms[subdomain_ID] = std::move(ho_ptr);
            context_vector[part]->hgo_parms[subdomain_ID] = std::move(hgo_ptr);
            context_vector[part]->neohookean_parms[subdomain_ID] = std::move(neohookean_ptr);
            context_vector[part]->exp_neohookean_parms[subdomain_ID] = std::move(exp_neohookean_ptr);
            context_vector[part]->biaxial_fung_type_parms[subdomain_ID] = std::move(biaxial_fung_type_ptr);
            context_vector[part]->augustin_parms[subdomain_ID] = std::move(augustin_ptr);
            context_vector[part]->nonlinear_spring_parms[subdomain_ID] = std::move(nonlinear_spring_ptr);
        }
    }

    // getting additional parameters for the heart part
    Pointer<Database> heart_parm_db = app_initializer->getComponentDatabase("ModelParameters");

    BoundaryConditions::t_load = heart_parm_db->getDouble("t_load");
    BoundaryConditions::P_load = heart_parm_db->getDouble("P_load");
    BoundaryConditions::kappa_surface_tether = heart_parm_db->getDouble("kappa_surface_tether");
    BoundaryConditions::kappa_body_tether = heart_parm_db->getDouble("kappa_body_tether");
    BoundaryConditions::kappa_surface_penalty = heart_parm_db->getDouble("kappa_surface_penalty");
    BoundaryConditions::kappa_body_penalty = heart_parm_db->getDouble("kappa_body_penalty");
    BoundaryConditions::kappa_pericardium = heart_parm_db->getDouble("kappa_pericardium");
    BoundaryConditions::eta_pericardium = heart_parm_db->getDouble("eta_pericardium");
    return;
}

void
ModelInitialization::create_part_dictionary()
{
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("AortaPart", &MeshInfo::AORTA));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("HeartPart", &MeshInfo::HEART));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("PulmArtPart", &MeshInfo::PULM_ART));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("AortaCapPart", &MeshInfo::AORTA_CAP));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("PulmArtCapPart", &MeshInfo::PULM_ART_CAP));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("AorticValvePart", &MeshInfo::AORTIC_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("PulmValvePart", &MeshInfo::PULM_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("VeinCapsPart", &MeshInfo::VEIN_CAPS));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("MitralValvePart", &MeshInfo::MITRAL_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("MitralPapillaryPart", &MeshInfo::MITRAL_PAP));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("MitralChordsPart", &MeshInfo::MITRAL_CHORDS));
    MeshInfo::part_dictionary.insert(
        std::pair<std::string, int*>("MitralStrutChordsPart", &MeshInfo::MITRAL_STRUT_CHORDS));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("TriValvePart", &MeshInfo::TRI_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("TriPapillaryPart", &MeshInfo::TRI_PAP));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("TriChordsPart", &MeshInfo::TRI_CHORDS));
    MeshInfo::part_dictionary.insert(
        std::pair<std::string, int*>("PulmArtWithValvePart", &MeshInfo::PULM_ART_WITH_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("AortaWithValvePart", &MeshInfo::AORTA_WITH_VALVE));
    MeshInfo::part_dictionary.insert(std::pair<std::string, int*>("EverythingElsePart", &MeshInfo::EVERYTHING_ELSE));
}
