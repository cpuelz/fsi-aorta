// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <tbox/NullDatabase.h>

// Headers for basic libMesh objects
#include <libmesh/analytic_function.h>
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/gmv_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/point.h>
#include <libmesh/quadrature_grid.h>
#include <libmesh/replicated_mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEInstrumentPanel.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFEPostProcessor.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// tools for initializing the finite element mesh
// and storing information about it
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/ModelInitialization.h>
#include <four_chambered_heart/PartContext.h>

// for the mechanics model
#include <four_chambered_heart/BoundaryConditions.h>
#include <four_chambered_heart/MechanicsModel.h>

// includes for afterload models
#include <four_chambered_heart/CirculationModel.h>
#include <four_chambered_heart/FeedbackForcer.h>
#include <four_chambered_heart/VelocityBcCoefs.h>

// includes for fluid sources
#include <four_chambered_heart/Source.h>
#include <four_chambered_heart/SourceDistributer.h>

// includes for model parameters
#include <four_chambered_heart/ModelParameters.h>

// C++ includes
#include <sstream>

// memory usage headers
#include <locale>
#include <sys/resource.h>
#include <sys/time.h>

// Fiber equation systems info
namespace FiberInfo
{
static const int FIBER_SYSTEM_ID = 0;
static const int FIBERX_ID = 0;
static const int FIBERY_ID = 1;
static const int FIBERZ_ID = 2;
} // namespace FiberInfo

using namespace FiberInfo;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 MeshBase& mesh,
                 EquationSystems* equation_systems,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

// to look at J, on the heart mesh.
void J_fcn(double& F,
           const libMesh::TensorValue<double>& FF,
           const libMesh::Point& x,
           const libMesh::Point& X,
           libMesh::Elem* elem,
           const std::vector<const std::vector<double>*>& system_var_data,
           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
           double data_time,
           void* ctx);

// to visualize the active tension function on the meshes...
void active_tension_fcn(double& F,
                        const libMesh::TensorValue<double>& FF,
                        const libMesh::Point& x,
                        const libMesh::Point& X,
                        libMesh::Elem* elem,
                        const std::vector<const std::vector<double>*>& system_var_data,
                        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
                        double data_time,
                        void* ctx);

namespace FCHTimer
{
struct Timer
{
    typedef std::chrono::high_resolution_clock::time_point timePoint_Type;
    typedef std::chrono::high_resolution_clock clock_Type;
    typedef std::chrono::seconds seconds;
    typedef std::chrono::duration<double> duration_Type;

    void reset()
    {
        M_run = false;
        M_elapsed = duration_Type::zero();
    }
    void start()
    {
        M_run = true;
        M_start = clock_Type::now();
    }
    void stop()
    {
        if (M_run)
        {
            M_run = false;
            M_end = clock_Type::now();
            M_elapsed += M_end - M_start;
        }
    }
    duration_Type elapsed()
    {
        if (!M_run)
            return M_elapsed;
        else
        {
            duration_Type elapsed = M_elapsed;
            auto end = clock_Type::now();
            elapsed += end - M_start;
            return elapsed;
        }
    }
    void print(std::ostream& out = std::cout)
    {
        out << "Timer: elapsed time = " << elapsed().count() << "  s." << std::endl;
    }
    duration_Type M_elapsed;
    timePoint_Type M_start;
    timePoint_Type M_end;
    bool M_run;
};

} // namespace Timer

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/

int
main(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        // use timer?
        const bool use_timer = input_db->getBoolWithDefault("use_timer", false);
        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string exodus_filename = app_initializer->getExodusIIFilename();
        const string viz_dump_dirname = app_initializer->getComponentDatabase("Main")->getString("viz_dump_dirname");
        const string viz_dump_file_extension = app_initializer->getComponentDatabase("Main")->getStringWithDefault("viz_dump_file_extension", ".e");
        const string viz_output_prename = input_db->getString("viz_output_prename");

        const bool from_restart = RestartManager::getManager()->isFromRestart();
        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();
        
        // whether to log invariant data
        const bool log_invariant_data = input_db->getBoolWithDefault("log_invariant_data", false);

        // whether to log active tension data
        const bool log_active_tension = input_db->getBoolWithDefault("log_active_tension", false);
        
        // whether to use volumetric energy
        const bool do_pressure_stabilization = input_db->getBoolWithDefault("DO_PRESSURE_STABILIZATION", false);

        // whether to split PK1 stress into dil and dev components
        const bool split_PK1_stress = input_db->getBoolWithDefault("SPLIT_PK1_STRESS", false);

        // bring in the meshes corresponding to each part
        const unsigned int NUM_PARTS = input_db->getInteger("NUM_PARTS");
        std::vector<std::unique_ptr<libMesh::MeshBase> > mesh_vector;
        std::vector<std::unique_ptr<ExodusII_IO> > exodus_mesh_readers;
        std::vector<std::string> part_names;
        // create part dictionary
        ModelInitialization::create_part_dictionary();
        std::map<std::string, int*>::iterator pp;
        int part = 0;

        // looping over possible parts
        for (pp = MeshInfo::part_dictionary.begin(); pp != MeshInfo::part_dictionary.end(); ++pp)
        {
            // checking to see if part name is in input file.
            if (!app_initializer->getInputDatabase()->keyExists(pp->first)) continue;
            Pointer<Database> part_db = app_initializer->getComponentDatabase(pp->first);
            mesh_vector.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
            exodus_mesh_readers.emplace_back(new ExodusII_IO(*mesh_vector[part]));
            exodus_mesh_readers[part]->read(part_db->getString("MESH_FILENAME"));
            *MeshInfo::part_dictionary[pp->first] = part;
            part_names.push_back(pp->first);
            ++part;
        }

        std::vector<std::unique_ptr<PartContext> > context_vector;
        for (unsigned int i = 0; i < NUM_PARTS; ++i)
        {
            context_vector.emplace_back(new PartContext(i));
        }

        // some functions need a vector of plain MeshBase pointers:
        std::vector<libMesh::MeshBase*> mesh_ptrs;
        for (const auto& mesh_ptr : mesh_vector)
        {
            mesh_ptrs.emplace_back(mesh_ptr.get());
        }

        if (NUM_PARTS != part_names.size())
        {
            TBOX_ERROR(
                "NUM_PARTS = " << NUM_PARTS
                               << " does not equal number of parts "
                                  "(supposedly "
                               << part_names.size()
                               << ") listed in input file. Check to make sure"
                                  "the part database names in the input file correspond to the names in the model "
                                  "dictionary.");
        }

        // prepare meshes for use
        for (unsigned int part = 0; part < mesh_ptrs.size(); ++part)
        {
            mesh_ptrs[part]->prepare_for_use();
        }

        // get mesh info
        ModelInitialization::load_mesh_info(mesh_ptrs, app_initializer, part_names);

        // bring in model parameters
        ModelInitialization::load_model_parameters(app_initializer, part_names, context_vector);

        // transform and scale heart mesh
        ModelInitialization::transform_meshes(mesh_ptrs, app_initializer, part_names);

        // determine order of the elements in each mesh part
        pout << "\n";
        for (int part = 0; part < input_db->getInteger("NUM_PARTS"); ++part)
        {
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            if (strcasecmp(elem_order.c_str(), "FIRST") == 0)
            {
                // There is nothing we can really do if the input mesh for some crazy
                // reason has second-order elements so just leave things alone
                pout << "mesh for " << part_names[part] << " is using FIRST order elements"
                     << "\n";
            }
            if (strcasecmp(elem_order.c_str(), "SECOND") == 0)
            {
                mesh_vector[part]->all_second_order();
                pout << "mesh for " << part_names[part] << " is using SECOND order elements"
                     << "\n";
            }
        }
        pout << "\n";
        
        if (from_restart)
        {
            // optionally uniformly refine meshes
            ModelInitialization::uniformly_refine_meshes(mesh_ptrs, app_initializer, part_names);
        
            // optionally locally refine meshes
            ModelInitialization::locally_refine_meshes(mesh_ptrs, app_initializer, part_names);
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }

        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           mesh_ptrs,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // create fiber_sys_data object for fiber fields
        pout << "**************************************************************************** \n";
        pout << "Setting up fiber info for PK1 stress...." << std::endl;
        pout << "**************************************************************************** \n";
        // this build a vector of SystemData of size 3 with NDIM=3
        // for the fiber and xfiber vectors.
        std::vector<int> vars1(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars1[d] = d;
        std::vector<SystemData> fiber_sys_data1(1);
        fiber_sys_data1[0] = SystemData(input_db->getString("FIBER_SYSTEM_NAME"), vars1);

        // this build a vector of SystemData of size 6 with NDIM=3
        // for the fiber and xfiber vectors.
        std::vector<int> vars2(2 * NDIM);
        for (unsigned int d = 0; d < 2 * NDIM; ++d) vars2[d] = d;
        std::vector<SystemData> fiber_sys_data2(1);
        fiber_sys_data2[0] = SystemData(input_db->getString("FIBER_SYSTEM_NAME"), vars2);

        // this build a vector of SystemData of size 9 with NDIM=3
        // for the fiber, sheets, and xfiber vectors.
        std::vector<int> vars3(3 * NDIM);
        for (unsigned int d = 0; d < 3 * NDIM; ++d) vars3[d] = d;
        std::vector<SystemData> fiber_sys_data3(1);
        fiber_sys_data3[0] = SystemData(input_db->getString("FIBER_SYSTEM_NAME"), vars3);

        // build system data vector for velocity
        std::vector<int> varv(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) varv[d] = d;
        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, varv);

        // Configure the IBFE solver.
        const std::string PK1_quad_order = input_db->getString("PK1_QUAD_ORDER");

        // things for everything else
        if (MeshInfo::EVERYTHING_ELSE >= 0)
        {
            const int part = MeshInfo::EVERYTHING_ELSE;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_everything_else(
                BoundaryConditions::tether_surface_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_everything_else(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_everything_else(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_everything_else.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_everything_else.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_everything_else, MeshInfo::EVERYTHING_ELSE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_everything_else, MeshInfo::EVERYTHING_ELSE);
            ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data_everything_else, MeshInfo::EVERYTHING_ELSE);
        }

        // things for the heart
        if (MeshInfo::HEART >= 0)
        {
            const int part = MeshInfo::HEART;
            IBFEMethod::LagBodyForceFcnData body_fcn_data_heart(BoundaryConditions::tether_body_force_function,
                                                                std::vector<IBTK::SystemData>(),
                                                                context_vector[part].get());
            IBFEMethod::LagSurfacePressureFcnData surface_pressure_data_heart(
                BoundaryConditions::loading_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::LagSurfaceForceFcnData surface_force_data_heart(
                BoundaryConditions::tether_surface_force_function, velocity_data, context_vector[part].get());
            ib_method_ops->registerLagSurfaceForceFunction(surface_force_data_heart, MeshInfo::HEART);
            ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data_heart, MeshInfo::HEART);
            ib_method_ops->registerLagBodyForceFunction(body_fcn_data_heart, MeshInfo::HEART);
            if (!do_pressure_stabilization)
            {
                if(split_PK1_stress)
                {
                    IBFEMethod::PK1StressFcnData PK1_dev_stress_data_heart(
                        MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    IBFEMethod::PK1StressFcnData PK1_dil_stress_data_heart(
                        MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    PK1_dev_stress_data_heart.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                    PK1_dil_stress_data_heart.quad_order = 
                        Utility::string_to_enum<libMesh::Order>(input_db->getString("DIL_QUAD_ORDER"));
                    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_heart, MeshInfo::HEART);
                    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_heart, MeshInfo::HEART);
                }
                else
                {
                    IBFEMethod::PK1StressFcnData PK1_combined_stress_data_heart(
                        MechanicsModel::PK1_combined_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    PK1_combined_stress_data_heart.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                    ib_method_ops->registerPK1StressFunction(PK1_combined_stress_data_heart, MeshInfo::HEART);
                }
            }
            else
            {
                IBFEMethod::PK1StressFcnData PK1_dev_stress_data_heart(
                    MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                PK1_dev_stress_data_heart.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_heart, MeshInfo::HEART);
            }
        }

        // things for the aorta
        if (MeshInfo::AORTA >= 0)
        {
             const int part = MeshInfo::AORTA;
            IBFEMethod::LagBodyForceFcnData body_fcn_data_aorta(BoundaryConditions::tether_body_force_function,
                                                                std::vector<IBTK::SystemData>(),
                                                                context_vector[part].get());
            IBFEMethod::LagSurfacePressureFcnData surface_pressure_data_aorta(
                BoundaryConditions::loading_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::LagSurfaceForceFcnData surface_force_data_aorta(
                BoundaryConditions::tether_surface_force_function, velocity_data, context_vector[part].get());
            ib_method_ops->registerLagSurfaceForceFunction(surface_force_data_aorta, MeshInfo::AORTA);
            ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data_aorta, MeshInfo::AORTA);
            ib_method_ops->registerLagBodyForceFunction(body_fcn_data_aorta, MeshInfo::AORTA);
            if (!do_pressure_stabilization)
            {
                if(split_PK1_stress)
                {
                    IBFEMethod::PK1StressFcnData PK1_dev_stress_data_aorta(
                        MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    IBFEMethod::PK1StressFcnData PK1_dil_stress_data_aorta(
                        MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    PK1_dev_stress_data_aorta.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                    PK1_dil_stress_data_aorta.quad_order = 
                        Utility::string_to_enum<libMesh::Order>(input_db->getString("DIL_QUAD_ORDER"));
                    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_aorta, MeshInfo::AORTA);
                    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_aorta, MeshInfo::AORTA);
                }
                else
                {
                    IBFEMethod::PK1StressFcnData PK1_combined_stress_data_aorta(
                        MechanicsModel::PK1_combined_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                    PK1_combined_stress_data_aorta.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                    ib_method_ops->registerPK1StressFunction(PK1_combined_stress_data_aorta, MeshInfo::AORTA);
                }
            }
            else
            {
                IBFEMethod::PK1StressFcnData PK1_dev_stress_data_aorta(
                    MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
                PK1_dev_stress_data_aorta.quad_order = Utility::string_to_enum<libMesh::Order>(PK1_quad_order);
                ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_aorta, MeshInfo::AORTA);
            }
        }

        // things for the PA
        if (MeshInfo::PULM_ART >= 0)
        {
            const int part = MeshInfo::PULM_ART;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_pa(BoundaryConditions::tether_surface_force_function,
                                                                   std::vector<IBTK::SystemData>(),
                                                                   context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_pa(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_pa(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_pa.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_pa.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_pa, MeshInfo::PULM_ART);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_pa, MeshInfo::PULM_ART);
            ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data_pa, MeshInfo::PULM_ART);
        }

        // things for the aortic valve
        if (MeshInfo::AORTIC_VALVE >= 0)
        {
            const int part = MeshInfo::AORTIC_VALVE;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_aortic_valve(
                BoundaryConditions::tether_surface_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_aortic_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data2, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_aortic_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_aortic_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_aortic_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_aortic_valve, MeshInfo::AORTIC_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_aortic_valve, MeshInfo::AORTIC_VALVE);
        }

        // things for the pulmonary valve
        if (MeshInfo::PULM_VALVE >= 0)
        {
            const int part = MeshInfo::PULM_VALVE;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_pulm_valve(
                BoundaryConditions::tether_surface_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_pulm_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data2, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_pulm_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_pulm_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_pulm_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_pulm_valve, MeshInfo::PULM_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_pulm_valve, MeshInfo::PULM_VALVE);
            // ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data_pulm_valve, MeshInfo::PULM_VALVE);
        }

        // things for the vein caps
        if (MeshInfo::VEIN_CAPS >= 0)
        {
            const int part = MeshInfo::VEIN_CAPS;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_vein_caps(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_vein_caps(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_vein_caps.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_vein_caps.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_vein_caps, MeshInfo::VEIN_CAPS);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_vein_caps, MeshInfo::VEIN_CAPS);
        }

        // things for the aorta cap
        if (MeshInfo::AORTA_CAP >= 0)
        {
            const int part = MeshInfo::AORTA_CAP;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_aorta_cap(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_aorta_cap(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_aorta_cap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_aorta_cap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_aorta_cap, MeshInfo::AORTA_CAP);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_aorta_cap, MeshInfo::AORTA_CAP);
        }

        // things for the pulm art cap
        if (MeshInfo::PULM_ART_CAP >= 0)
        {
            const int part = MeshInfo::PULM_ART_CAP;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_pulm_art_cap(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_pulm_art_cap(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_pulm_art_cap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_pulm_art_cap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_pulm_art_cap, MeshInfo::PULM_ART_CAP);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_pulm_art_cap, MeshInfo::PULM_ART_CAP);
        }

        // things for the mitral valve
        if (MeshInfo::MITRAL_VALVE >= 0)
        {
            const int part = MeshInfo::MITRAL_VALVE;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_mitral_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data2, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_mitral_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_mitral_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_mitral_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_mitral_valve, MeshInfo::MITRAL_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_mitral_valve, MeshInfo::MITRAL_VALVE);
        }

        // things for the mitral papillary muscles
        if (MeshInfo::MITRAL_PAP >= 0)
        {
            const int part = MeshInfo::MITRAL_PAP;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_mitral_pap(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data1, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_mitral_pap(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_mitral_pap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_mitral_pap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_mitral_pap, MeshInfo::MITRAL_PAP);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_mitral_pap, MeshInfo::MITRAL_PAP);
        }

        // things for the mitral chordae
        if (MeshInfo::MITRAL_CHORDS >= 0)
        {
            const int part = MeshInfo::MITRAL_CHORDS;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_mitral_chords(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data1, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_mitral_chords(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_mitral_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_mitral_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_mitral_chords, MeshInfo::MITRAL_CHORDS);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_mitral_chords, MeshInfo::MITRAL_CHORDS);
        }

        // things for the mitral strut chordae
        if (MeshInfo::MITRAL_STRUT_CHORDS >= 0)
        {
            const int part = MeshInfo::MITRAL_STRUT_CHORDS;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_mitral_strut_chords(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data1, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_mitral_strut_chords(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_mitral_strut_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_mitral_strut_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_mitral_strut_chords,
                                                     MeshInfo::MITRAL_STRUT_CHORDS);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_mitral_strut_chords,
                                                     MeshInfo::MITRAL_STRUT_CHORDS);
        }

        // things for the tricuspid valve
        if (MeshInfo::TRI_VALVE >= 0)
        {
            const int part = MeshInfo::TRI_VALVE;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_tri_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data2, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_tri_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_tri_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_tri_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_tri_valve, MeshInfo::TRI_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_tri_valve, MeshInfo::TRI_VALVE);
        }

        // things for the tricuspid papillary muscles
        if (MeshInfo::TRI_PAP >= 0)
        {
            const int part = MeshInfo::TRI_PAP;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_tri_pap(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data1, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_tri_pap(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_tri_pap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_tri_pap.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_tri_pap, MeshInfo::TRI_PAP);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_tri_pap, MeshInfo::TRI_PAP);
        }

        // things for the tricuspid chordae
        if (MeshInfo::TRI_CHORDS >= 0)
        {
            const int part = MeshInfo::TRI_CHORDS;
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_tri_chords(
                MechanicsModel::PK1_dev_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_tri_chords(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_tri_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_tri_chords.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_tri_chords, MeshInfo::TRI_CHORDS);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_tri_chords, MeshInfo::TRI_CHORDS);
        }

        // things for the aorta with valve
        if (MeshInfo::AORTA_WITH_VALVE >= 0)
        {
            const int part = MeshInfo::AORTA_WITH_VALVE;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_aorta_with_valve(
                BoundaryConditions::tether_surface_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_aorta_with_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data2, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_aorta_with_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_aorta_with_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_aorta_with_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_aorta_with_valve, MeshInfo::AORTA_WITH_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_aorta_with_valve, MeshInfo::AORTA_WITH_VALVE);
            ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data_aorta_with_valve,
                                                           MeshInfo::AORTA_WITH_VALVE);
        }

        // things for the pulmonary artery with valve
        if (MeshInfo::PULM_ART_WITH_VALVE >= 0)
        {
            const int part = MeshInfo::PULM_ART_WITH_VALVE;
            IBFEMethod::LagSurfaceForceFcnData surface_fcn_data_pulm_art_with_valve(
                BoundaryConditions::tether_surface_force_function,
                std::vector<IBTK::SystemData>(),
                context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dev_stress_data_pulm_art_with_valve(
                MechanicsModel::PK1_dev_stress_function, fiber_sys_data3, context_vector[part].get());
            IBFEMethod::PK1StressFcnData PK1_dil_stress_data_pulm_art_with_valve(
                MechanicsModel::PK1_dil_stress_function, std::vector<IBTK::SystemData>(), context_vector[part].get());
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            const std::string elem_order = part_db->getString("ELEM_ORDER");
            PK1_dev_stress_data_pulm_art_with_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            PK1_dil_stress_data_pulm_art_with_valve.quad_order = Utility::string_to_enum<libMesh::Order>(
                strcasecmp(elem_order.c_str(), "SECOND") == 0 ? "FIFTH" : "THIRD");
            ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_pulm_art_with_valve,
                                                     MeshInfo::PULM_ART_WITH_VALVE);
            ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_pulm_art_with_valve,
                                                     MeshInfo::PULM_ART_WITH_VALVE);
            ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data_pulm_art_with_valve,
                                                           MeshInfo::PULM_ART_WITH_VALVE);
        }

        // initialize libmesh systems
        ib_method_ops->initializeFEEquationSystems();
        
        // pressure stabilization
        if(do_pressure_stabilization)
        {
        //    ib_method_ops->registerPressureStabilizationPart(MeshInfo::HEART);
        }

        // populate equation system vector for each part
        std::vector<EquationSystems*> eq_systems_vector;
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            eq_systems_vector.push_back(ib_method_ops->getFEDataManager(part)->getEquationSystems());
        }

        // initializing post processors to reconstruct J
        std::vector<std::unique_ptr<IBFECentroidPostProcessor> > ib_post_processors;
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            pout << "Configuring the postprocessor...\n";
            ib_post_processors.emplace_back(
                new IBFECentroidPostProcessor("IBFEPostProcessor", ib_method_ops->getFEDataManager(part)));
            // set up evaluation of the determinant of the deformation gradient
            ib_post_processors[part]->registerScalarVariable(
                "J", MONOMIAL, CONSTANT, &J_fcn, std::vector<IBTK::SystemData>(), context_vector[part].get());
            ib_post_processors[part]->registerScalarVariable("active tension",
                                                             MONOMIAL,
                                                             CONSTANT,
                                                             &active_tension_fcn,
                                                             std::vector<IBTK::SystemData>(),
                                                             context_vector[part].get());
        }

        // Create Eulerian boundary condition specification objects.
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        std::unique_ptr<CirculationModel> circ_model;
        if (MeshInfo::afterload_parms.size() > 0)
        // if this vector is not zero length, then afterload models exist, so set them up.
        {
            pout << "\n**************************************************************************** \n";
            pout << "setting up afterload models by constructing CirculationModel object\n";
            pout << "**************************************************************************** \n";
            circ_model.reset(new CirculationModel(
                "circ_model", input_db->getDatabase("AfterloadModel"), MeshInfo::afterload_parms, true));
            for (int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = new VelocityBcCoefs(circ_model.get(), d);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
            Pointer<FeedbackForcer> feedback_forcer =
                new FeedbackForcer(circ_model.get(), navier_stokes_integrator, patch_hierarchy);
            time_integrator->registerBodyForceFunction(feedback_forcer);
        }
        else // do something simpler
        {
            const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
            if (periodic_shift.min() <= 0)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    ostringstream bc_coefs_name_stream;
                    bc_coefs_name_stream << "u_bc_coefs_" << d;
                    const string bc_coefs_name = bc_coefs_name_stream.str();

                    ostringstream bc_coefs_db_name_stream;
                    bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                    const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                    u_bc_coefs[d] = new muParserRobinBcCoefs(
                        bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
            }
        }

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // set up exodus II IO objects for each part
        std::vector<std::unique_ptr<ExodusII_IO> > exodus_io_vector;
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            exodus_io_vector.emplace_back(uses_exodus ? new ExodusII_IO(*mesh_vector[part]) : nullptr);
        }

        // check to see if this is a restarted run.  if so, we just need to append to our
        // current exodus files
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            if (exodus_io_vector[part].get())
            {
                exodus_io_vector[part]->append(from_restart);
            }
        }

        // Initialize fiber systems (but do not load the data yet).
        if (!from_restart)
        {
            for (unsigned int part = 0; part < mesh_vector.size(); ++part)
            {
                if (MeshInfo::has_fibers[part])
                {
                    pout << "\n";
                    pout << "**************************************************************************** \n";
                    pout << "Initializing fiber systems for mesh part " << part << ": " << part_names[part] << "\n";
                    ModelInitialization::init_fiber_system(
                    *exodus_mesh_readers[part],
                            *eq_systems_vector[part],
                            input_db,
                            ib_method_ops->getPK1StressFunction(part)[0].system_data[0].vars.size());
                    pout << "**************************************************************************** \n";
                    pout << "\n";
                }
            }
        }

        // NOTE: All libMesh Systems should be set up before this point in the routine.
        ib_method_ops->initializeFEData();
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            if (ib_post_processors[part]) ib_post_processors[part]->initializeFEData();
        }

        // initialize fluid sources
        std::vector<std::unique_ptr<Source> > source_vector;
        for (unsigned int part = 0; part < mesh_vector.size(); ++part)
        {
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            if (MeshInfo::turn_on_sources[part])
            {
                pout << "\n";
                pout << "**************************************************************************** \n";
                pout << "Initializing fluid sources for mesh part " << part << ": " << part_names[part] << "\n";
                SAMRAI::tbox::Array<std::string> source_names = part_db->getStringArray("source_names");
                for (int ss = 0; ss < source_names.size(); ++ss)
                {
                    const std::string source_db_name = "source_" + source_names[ss];
                    Pointer<Database> source_db = app_initializer->getComponentDatabase(source_db_name);
                    source_vector.emplace_back(
                        new Source(source_names[ss], part, source_db, ib_method_ops, true /*register_for_restart*/));
                }
                pout << "**************************************************************************** \n";
                pout << "\n";
            }
        }
        Pointer<SourceDistributer> source_distributor =
            (source_vector.size() > 0) ?
                new SourceDistributer(std::move(source_vector), navier_stokes_integrator, patch_hierarchy) :
                NULL;
        if (source_distributor) navier_stokes_integrator->registerFluidSourceFunction(source_distributor);

        // Initialize hierarchy configuration and data on all patches.
        //
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Read in the fiber data.
        for (unsigned int part = 0; part < mesh_vector.size(); ++part)
        {
            if (MeshInfo::has_fibers[part])
            {
                pout << "\n";
                pout << "**************************************************************************** \n";
                pout << "Copying fiber information from file into equation systems object for mesh part " << part
                     << ": " << part_names[part] << "\n";
                if (!from_restart)
                {
                    ModelInitialization::load_fiber_system(*exodus_mesh_readers[part], *eq_systems_vector[part], input_db);
                }
                pout << "**************************************************************************** \n";
                pout << "\n";
            }
        }

        if (!from_restart)
        {
            // optionally uniformly refine meshes
            ModelInitialization::uniformly_refine_meshes(mesh_ptrs, app_initializer, part_names);
        
            // optionally locally refine meshes
            ModelInitialization::locally_refine_meshes(mesh_ptrs, app_initializer, part_names);
            
            // reinit equation systems
            ModelInitialization::reinit_equation_systems(eq_systems_vector, app_initializer, part_names);
        }

        // We have a bit of a catch-22 when it comes to mesh refinement:
        // 1. IBFEMethod assumes that the mesh does not change after the call to
        //    initializeFEData().
        // 2. We cannot load the fiber fields until the solution vectors exist
        //    (which happens after dofs are distributed in initializeFEData).
        // 3. We have to load the fiber fields before doing mesh refinement so
        //    that they are prolongated correctly.
        //
        // Hence we have to do things that are not yet officially supported in
        // IBFEMethod - i.e., reinitialize all the finite element data at this
        // point, to get around the first problem.
        //
        // This fix works because we only need to read some data from the
        // FEDataManager objects in the postprocessor below: we can rely on
        // the initial regrid to set up everything else we will use later.
        for (unsigned int part = 0; part < mesh_vector.size(); ++part)
        {
            ib_method_ops->getFEDataManager(part)->reinitElementMappings();
        }
        
        // initialize IBFE instrumentation
        std::vector<std::unique_ptr<IBFEInstrumentPanel> > instrument_vector;
        std::map<int, int> part_ID_to_inst_ID;
        int ii = 0;
        for (unsigned int part = 0; part < mesh_vector.size(); ++part)
        {
            Pointer<Database> part_db = app_initializer->getComponentDatabase(part_names[part]);
            if (MeshInfo::turn_on_meters[part])
            {
                instrument_vector.emplace_back(new IBFEInstrumentPanel(part_db, part));
                part_ID_to_inst_ID.insert(std::pair<int, int>(part, ii));
                ii += 1;
            }
        }

        // initialize data for meters
        for (unsigned int part = 0; part < mesh_vector.size(); ++part)
        {
            if (MeshInfo::turn_on_meters[part])
            {
                instrument_vector[part_ID_to_inst_ID[part]]->initializeHierarchyIndependentData(ib_method_ops);
            }
        }

        // Deallocate initialization objects.
        app_initializer.setNull();

        // clear the exodus readers.
        exodus_mesh_readers.clear();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                for (unsigned int part = 0; part < NUM_PARTS; ++part)
                {
                    if (ib_post_processors[part]) ib_post_processors[part]->postProcessData(loop_time);
                    exodus_io_vector[part]->write_timestep(viz_dump_dirname + "/" + part_names[part] + viz_dump_file_extension,
                                                           *eq_systems_vector[part],
                                                           iteration_num / viz_dump_interval + 1,
                                                           loop_time);
                    if (MeshInfo::turn_on_meters[part])
                    {
                        const std::string meter_plot_dir_name =
                            instrument_vector[part_ID_to_inst_ID[part]]->getPlotDirectoryName();
                        const unsigned int num_meters =
                            instrument_vector[part_ID_to_inst_ID[part]]->getNumberOfMeterMeshes();
                        for (unsigned int jj = 0; jj < num_meters; ++jj)
                        {
                            const std::string meter_mesh_name =
                                instrument_vector[part_ID_to_inst_ID[part]]->getMeterMeshName(jj);
                            libMesh::ExodusII_IO meter_exodus_io(
                                instrument_vector[part_ID_to_inst_ID[part]]->getMeterMesh(jj));
                            libMesh::EquationSystems& meter_system =
                                instrument_vector[part_ID_to_inst_ID[part]]->getMeterMeshEquationSystems(jj);
                            std::ostringstream mesh_output;
                            mesh_output
                                << meter_plot_dir_name
                                << "/"
                                << "" << meter_mesh_name << viz_dump_file_extension;
                            meter_exodus_io.write_timestep(
                                mesh_output.str(), meter_system, iteration_num / viz_dump_interval + 1, loop_time);
                        }
                    }
                }
            }
        }

        // Open streams to save volume of structure.
        ofstream volume_stream;
        ofstream J_stream;
        ofstream I1_stream;
        const std::string file_prename = input_db->getString("file_prename");
        if (SAMRAI_MPI::getRank() == 0)
        {
            std::ostringstream volume_output;
            std::ostringstream J_output;
            std::ostringstream I1_output;
            volume_output << viz_output_prename + "volume_curve.dat";
            J_output << viz_output_prename + "J_values.dat";
            I1_output << viz_output_prename + "I1_values.dat";
            volume_stream.open(volume_output.str().c_str(), ios_base::out | ios_base::trunc);
            J_stream.open(J_output.str().c_str(), ios_base::out | ios_base::trunc);
            I1_stream.open(I1_output.str().c_str(), ios_base::out | ios_base::trunc);
        }

        //********************************************************************************
        // setting up some objects for measuring fluxes and mean pressures
        //********************************************************************************
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
        Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        Pointer<VariableContext> u_current_ctx = navier_stokes_integrator->getCurrentContext();
        const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, p_current_ctx);
        const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, u_current_ctx);
        const double current_time = 0.0;
        Pointer<SideVariable<NDIM, double> > u_copy_var = new SideVariable<NDIM, double>("u_copy");
        Pointer<CellVariable<NDIM, double> > p_copy_var = new CellVariable<NDIM, double>("p_copy");
        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const int u_copy_idx =
            var_db->registerVariableAndContext(u_copy_var, time_integrator->getScratchContext(), ib_ghosts);
        const int p_copy_idx =
            var_db->registerVariableAndContext(p_copy_var, time_integrator->getScratchContext(), ib_ghosts);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_copy_idx, current_time);
            level->allocatePatchData(p_copy_idx, current_time);
        }

        // **********************************************
        // get mean pressure and velocity on surface mesh
        //***********************************************
        HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
            hier_data_ops_manager->getOperationsDouble(p_var, patch_hierarchy, true);
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_sc_data_ops =
            hier_data_ops_manager->getOperationsDouble(u_var, patch_hierarchy, true);
        hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);
        hier_sc_data_ops->copyData(u_copy_idx, u_current_idx, true);

        // fill ghost cells for pressure
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> p_transaction_comp(1);
        p_transaction_comp[0] = InterpolationTransactionComponent(p_copy_idx,
                                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                                  /*use_cf_bdry_interpolation*/ false,
                                                                  "CONSERVATIVE_COARSEN",
                                                                  "LINEAR");

        Pointer<HierarchyGhostCellInterpolation> p_hier_bdry_fill = new HierarchyGhostCellInterpolation();
        p_hier_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
        p_hier_bdry_fill->fillData(loop_time);

        // fill ghost cells for velocity
        std::vector<InterpolationTransactionComponent> u_transaction_comp(1);
        u_transaction_comp[0] = InterpolationTransactionComponent(u_copy_idx,
                                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                                  /*use_cf_bdry_interpolation*/ false,
                                                                  "CONSERVATIVE_COARSEN",
                                                                  "LINEAR");

        Pointer<HierarchyGhostCellInterpolation> u_hier_bdry_fill = new HierarchyGhostCellInterpolation();
        u_hier_bdry_fill->initializeOperatorState(u_transaction_comp, patch_hierarchy);
        u_hier_bdry_fill->fillData(loop_time);

        // read instrument data
        for (unsigned int part = 0; part < NUM_PARTS; ++part)
        {
            if (MeshInfo::turn_on_meters[part])
            {
                instrument_vector[part_ID_to_inst_ID[part]]->initializeHierarchyDependentData(ib_method_ops,
                                                                                              patch_hierarchy);
                instrument_vector[part_ID_to_inst_ID[part]]->readInstrumentData(
                    u_copy_idx, p_copy_idx, patch_hierarchy, loop_time);
            }
        }
        //********************************************************************************/

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        int count_output = 0;
        bool inverted_element_check = false;

        // Create Global Timer

        FCHTimer::Timer timeloop_timer;
        FCHTimer::Timer timestep_timer;
        FCHTimer::Timer check_J_timer;

        if (use_timer)
        {
            SAMRAI_MPI::barrier();
            timeloop_timer.start();
        }

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            // for checking the circ models
            if (circ_model->d_test_circ_models)
            {
                dt = time_integrator->getMaximumTimeStepSize();
                loop_time += dt;
                std::cout << "time = " << loop_time << "\n";
                Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                Pointer<HierarchyMathOps> hier_math_ops_foo = navier_stokes_integrator->getHierarchyMathOps();
                const int wgt_cc_idx = hier_math_ops_foo->getCellWeightPatchDescriptorIndex();
                const int wgt_sc_idx = hier_math_ops_foo->getSideWeightPatchDescriptorIndex();
                if (circ_model)
                {
                    circ_model->advanceTimeDependentData(
                    dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
                }
            }
            // for solving the entire model!
            else
            {
                rusage ru;
                getrusage(RUSAGE_SELF, &ru);
                long mem = ru.ru_maxrss;
                pout << "Current memory max:" << long(SAMRAI_MPI::sumReduction(double(mem))) << "\n";
                
                if (use_timer)
                {
                    timestep_timer.reset();
                    SAMRAI_MPI::barrier();
                    timestep_timer.start();
                }
                
                iteration_num = time_integrator->getIntegratorStep();
                loop_time = time_integrator->getIntegratorTime();
                
                // check to see if we have an inverted element (or one close to inverting)
                inverted_element_check = SAMRAI_MPI::maxReduction(MechanicsModel::inverted_element ? 1 : 0) > 0;
                
                // Compute min and max of J, I1, and I4's for each subdomain of each part
                for (unsigned int part = 0; part < NUM_PARTS; ++part)
                {
                    for (unsigned int ss = 0; ss < MeshInfo::subdomain_IDs[part].size(); ++ss)
                    {
                        const libMesh::subdomain_id_type subdomain_ID = MeshInfo::subdomain_IDs[part][ss];
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
                        context_vector[part]->active_tension[subdomain_ID] = -0.5 * std::numeric_limits<double>::max();
                    }
                }
                
                pout << "\n";
                pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                pout << "At beginning of timestep # " << iteration_num << "\n";
                pout << "Simulation time is " << loop_time << "\n";
                
                dt = time_integrator->getMaximumTimeStepSize();
                
                Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                Pointer<HierarchyMathOps> hier_math_ops_foo = navier_stokes_integrator->getHierarchyMathOps();
                const int wgt_cc_idx = hier_math_ops_foo->getCellWeightPatchDescriptorIndex();
                const int wgt_sc_idx = hier_math_ops_foo->getSideWeightPatchDescriptorIndex();
                if (circ_model)
                {
                    circ_model->advanceTimeDependentData(
                    dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
                }
                
                time_integrator->advanceHierarchy(dt);
                loop_time += dt;
                
                // update source pressures.
                if (source_distributor)
                {
                    source_distributor->interpolatePressure(P_current_idx, loop_time);
                    source_distributor->updateSources(ib_method_ops, loop_time, dt);
                }
                
                pout << "\n";
                pout << "At end       of timestep # " << iteration_num << "\n";
                pout << "Simulation time is " << loop_time << "\n";
                pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                pout << "\n";
                
                // see if we need to reallocate
                HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
                hier_math_ops.setPatchHierarchy(patch_hierarchy);
                hier_math_ops.resetLevels(coarsest_ln, finest_ln);
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    if (!level->checkAllocated(p_copy_idx)) level->allocatePatchData(p_copy_idx);
                    if (!level->checkAllocated(u_copy_idx)) level->allocatePatchData(u_copy_idx);
                }
                
                hier_cc_data_ops->copyData(p_copy_idx, p_current_idx, true);
                hier_sc_data_ops->copyData(u_copy_idx, u_current_idx, true);
                
                p_hier_bdry_fill->initializeOperatorState(p_transaction_comp, patch_hierarchy);
                p_hier_bdry_fill->fillData(loop_time);
                
                u_hier_bdry_fill->initializeOperatorState(u_transaction_comp, patch_hierarchy);
                u_hier_bdry_fill->fillData(loop_time);
                
                if (use_timer)
                {
                    SAMRAI_MPI::barrier();
                    timestep_timer.stop();
                    timestep_timer.print(pout);
                }
                // write out data from meter meshes
                for (unsigned int part = 0; part < NUM_PARTS; ++part)
                {
                    if (MeshInfo::turn_on_meters[part])
                    {
                        if ((iteration_num % instrument_vector[part_ID_to_inst_ID[part]]->getInstrumentDumpInterval()) == 0)
                        {
                            instrument_vector[part_ID_to_inst_ID[part]]->initializeHierarchyDependentData(ib_method_ops,
                                    patch_hierarchy);
                            instrument_vector[part_ID_to_inst_ID[part]]->readInstrumentData(
                            u_copy_idx, p_copy_idx, patch_hierarchy, loop_time);
                        }
                    }
                }
                if (use_timer)
                {
                    check_J_timer.reset();
                    SAMRAI_MPI::barrier();
                    check_J_timer.start();
                }
                
                if (log_invariant_data)
                {
                    for (unsigned int part = 0; part < NUM_PARTS; ++part)
                    {
                        
                        plog << "-------------------\n";
                        plog << "Data for part " << part << ": \n";
                        for (unsigned int ss = 0; ss < MeshInfo::subdomain_IDs[part].size(); ++ss)
                        {
                            const libMesh::subdomain_id_type subdomain_ID = MeshInfo::subdomain_IDs[part][ss];
                            plog << "subdomain " << subdomain_ID << " with material model name ''"
                                    << context_vector[part]->constitutive_model_name[subdomain_ID] << "'':\n";
                            context_vector[part]->min_J[subdomain_ID] =
                                    SAMRAI_MPI::minReduction(context_vector[part]->min_J[subdomain_ID]);
                            context_vector[part]->max_J[subdomain_ID] =
                                    SAMRAI_MPI::maxReduction(context_vector[part]->max_J[subdomain_ID]);
                            context_vector[part]->min_I1[subdomain_ID] =
                                    SAMRAI_MPI::minReduction(context_vector[part]->min_I1[subdomain_ID]);
                            context_vector[part]->max_I1[subdomain_ID] =
                                    SAMRAI_MPI::maxReduction(context_vector[part]->max_I1[subdomain_ID]);
                            plog << "min/max J = " << context_vector[part]->min_J[subdomain_ID] << "/"
                                    << context_vector[part]->max_J[subdomain_ID] << "\n";
                            plog << "min/max I1 = " << context_vector[part]->min_I1[subdomain_ID] << "/"
                                    << context_vector[part]->max_I1[subdomain_ID] << "\n";
                            if (MeshInfo::has_fibers[part])
                            {
                                context_vector[part]->min_I4_one[subdomain_ID] =
                                        SAMRAI_MPI::minReduction(context_vector[part]->min_I4_one[subdomain_ID]);
                                context_vector[part]->max_I4_one[subdomain_ID] =
                                        SAMRAI_MPI::maxReduction(context_vector[part]->max_I4_one[subdomain_ID]);
                                context_vector[part]->min_I4_two[subdomain_ID] =
                                        SAMRAI_MPI::minReduction(context_vector[part]->min_I4_two[subdomain_ID]);
                                context_vector[part]->max_I4_two[subdomain_ID] =
                                        SAMRAI_MPI::maxReduction(context_vector[part]->max_I4_two[subdomain_ID]);
                                context_vector[part]->min_I4_three[subdomain_ID] =
                                        SAMRAI_MPI::minReduction(context_vector[part]->min_I4_three[subdomain_ID]);
                                context_vector[part]->max_I4_three[subdomain_ID] =
                                        SAMRAI_MPI::maxReduction(context_vector[part]->max_I4_three[subdomain_ID]);
                                plog << "min/max I4_one = " << context_vector[part]->min_I4_one[subdomain_ID] << "/"
                                        << context_vector[part]->max_I4_one[subdomain_ID] << "\n";
                                plog << "min/max I4_two = " << context_vector[part]->min_I4_two[subdomain_ID] << "/"
                                        << context_vector[part]->max_I4_two[subdomain_ID] << "\n";
                                plog << "min/max I4_three = " << context_vector[part]->min_I4_three[subdomain_ID] << "/"
                                        << context_vector[part]->max_I4_three[subdomain_ID] << "\n\n";
                            }
                            else
                            {
                                plog << "\n";
                            }
                        }
                        plog << "-------------------\n";
                    }
                }
               
                if (log_active_tension)
                {
                    for (unsigned int part = 0; part < NUM_PARTS; ++part)
                    {
                        
                        plog << "-------------------\n";
                        plog << "Active tension data for part " << part << ": \n";
                        for (unsigned int ss = 0; ss < MeshInfo::subdomain_IDs[part].size(); ++ss)
                        {
                            const libMesh::subdomain_id_type subdomain_ID = MeshInfo::subdomain_IDs[part][ss];
                            plog << "subdomain " << subdomain_ID << " with material model name ''"
                                    << context_vector[part]->constitutive_model_name[subdomain_ID] << "'':\n";
                            context_vector[part]->active_tension[subdomain_ID] =
                                    SAMRAI_MPI::maxReduction(context_vector[part]->active_tension[subdomain_ID]);
                            plog << "active tension = " << context_vector[part]->active_tension[subdomain_ID] << "\n";
                        }
                        plog << "-------------------\n";
                    }
                }
                
                // Compute min and max I1 and J for each part.
                if (inverted_element_check)
                {
                    for (unsigned int part = 0; part < NUM_PARTS; ++part)
                    {
                        System& X_system = eq_systems_vector[part]->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
                        NumericVector<double>* X_vec = X_system.solution.get();
                        NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
                        copy_and_synch(*X_vec, *X_ghost_vec);
                        double J_integral = 0.0;
                        DofMap& X_dof_map = X_system.get_dof_map();
                        std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
                        std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
                        std::unique_ptr<QBase> qrule;
                        if(split_PK1_stress)
                        {
                            qrule = QBase::build(QGAUSS,
                                                 NDIM, 
                                                 Utility::string_to_enum<libMesh::Order>(input_db->getString("DIL_QUAD_ORDER")));
                        }
                        else
                        {
                            qrule = QBase::build(QGAUSS, 
                                                 NDIM, 
                                                 Utility::string_to_enum<libMesh::Order>(PK1_quad_order));
                        }
                        fe->attach_quadrature_rule(qrule.get());
                        const std::vector<double>& JxW = fe->get_JxW();
                        const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
                        TensorValue<double> FF;
                        boost::multi_array<double, 2> X_node;
                        const MeshBase::const_element_iterator el_begin = mesh_vector[part]->active_local_elements_begin();
                        const MeshBase::const_element_iterator el_end = mesh_vector[part]->active_local_elements_end();
                        double max_J = -0.5 * std::numeric_limits<double>::max();
                        double min_J = 0.5 * std::numeric_limits<double>::max();
                        double max_I1 = -0.5 * std::numeric_limits<double>::max();
                        double min_I1 = 0.5 * std::numeric_limits<double>::max();
                        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
                        {
                            Elem* const elem = *el_it;
                            fe->reinit(elem);
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
                            }
                            const int n_qp = qrule->n_points();
                            get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
                            for (int qp = 0; qp < n_qp; ++qp)
                            {
                                jacobian(FF, qp, X_node, dphi);
                                J_integral += abs(FF.det()) * JxW[qp];
                                min_J = std::min(min_J, FF.det());
                                max_J = std::max(max_J, FF.det());
                                min_I1 = std::min(min_I1, FF.tr());
                                max_I1 = std::max(max_I1, FF.tr());
                            }
                        }
                        J_integral = SAMRAI_MPI::sumReduction(J_integral);
                        min_J = SAMRAI_MPI::minReduction(min_J);
                        max_J = SAMRAI_MPI::maxReduction(max_J);
                        min_I1 = SAMRAI_MPI::minReduction(min_I1);
                        max_I1 = SAMRAI_MPI::maxReduction(max_I1);
                        
                        if (SAMRAI_MPI::getRank() == 0)
                        {
                            volume_stream.precision(4);
                            volume_stream.setf(ios::scientific, ios::floatfield);
                            
                            J_stream.precision(3);
                            J_stream.setf(ios::scientific, ios::floatfield);
                            
                            I1_stream.precision(3);
                            I1_stream.setf(ios::scientific, ios::floatfield);
                            if (part == 0 && NUM_PARTS > 1)
                            {
                                I1_stream << loop_time << " " << min_I1 << " " << max_I1 << " ";
                                J_stream << loop_time << " " << min_J << " " << max_J << " ";
                            }
                            else if (part == NUM_PARTS - 1)
                            {
                                I1_stream << min_I1 << " " << max_I1 << endl;
                                J_stream << min_J << " " << max_J << endl;
                            }
                            else
                            {
                                I1_stream << min_I1 << " " << max_I1 << " ";
                                J_stream << min_J << " " << max_J << " ";
                            }
                        }
                    }
                } // if inverted element
                
                if (use_timer)
                {
                    SAMRAI_MPI::barrier();
                    check_J_timer.stop();
                    check_J_timer.print(pout);
                }
                // At specified intervals, write visualization and restart files,
                // print out timer data, and store hierarchy data for post
                // processing.
                iteration_num += 1;
                const bool last_step = !time_integrator->stepsRemaining();
                if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step || inverted_element_check))
                {
                    count_output += 1;
                    pout << "\nWriting visualization files...\n\n";
                    if (uses_visit)
                    {
                        time_integrator->setupPlotData();
                        visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                    }
                    if (uses_exodus)
                    {
                        for (unsigned int part = 0; part < NUM_PARTS; ++part)
                        {
                            if (ib_post_processors[part]) ib_post_processors[part]->postProcessData(loop_time);
                            exodus_io_vector[part]->write_timestep(viz_dump_dirname + "/" + part_names[part] + viz_dump_file_extension,
                                    *eq_systems_vector[part],
                                    iteration_num / viz_dump_interval + 1,
                                    loop_time);
                        }
                    }
                }
                if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
                {
                    pout << "\nWriting restart files...\n\n";
                    RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                    ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
                }
                if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
                {
                    pout << "\nWriting timer data...\n\n";
                    TimerManager::getManager()->print(plog);
                }
                if ((dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step)) &&
                        MeshInfo::HEART >= 0)
                {
                    pout << "\nWriting state data...\n\n";
                    output_data(patch_hierarchy,
                            navier_stokes_integrator,
                            *mesh_vector[MeshInfo::HEART],
                            eq_systems_vector[MeshInfo::HEART],
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
                }
                
                // check to see if we have an inverted element
                if (inverted_element_check)
                {
                    TBOX_ERROR(
                            "The determinant of the deformation gradient "
                            "at a quadrature point is close to zero or negative, "
                            "so an element may have inverted.");
                }
            }
        }

        if (use_timer)
        {
            SAMRAI_MPI::barrier();
            timeloop_timer.stop();
            timeloop_timer.print(pout);
        }
        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (auto& u_bc_coef : u_bc_coefs)
        {
            delete u_bc_coef;
        }
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}

void
J_fcn(double& F,
      const libMesh::TensorValue<double>& FF,
      const libMesh::Point& /*x*/,
      const libMesh::Point& /*X*/,
      libMesh::Elem* /*elem*/,
      const std::vector<const std::vector<double>*>& /*system_var_data*/,
      const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
      double /*data_time*/,
      void* /*ctx*/)
{
    F = FF.det();
    return;
}

void
active_tension_fcn(double& F,
                   const libMesh::TensorValue<double>& /*FF*/,
                   const libMesh::Point& /*x*/,
                   const libMesh::Point& /*X*/,
                   libMesh::Elem* elem,
                   const std::vector<const std::vector<double>*>& /*system_var_data*/,
                   const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
                   double data_time,
                   void* ctx)
{
    F = 0.0;
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const bool enable_active_stress = part_context->general_parms.at(subdomain_ID)->enable_active_stress;
    const bool enable_active_strain = part_context->general_parms.at(subdomain_ID)->enable_active_strain;

    if (enable_active_stress || enable_active_strain)
    {
        const double max_tension = part_context->general_parms.at(subdomain_ID)->Tension;
        F = max_tension * part_context->general_parms.at(subdomain_ID)->active_tension_function(data_time, elem, ctx);
    }
    return;
}

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            MeshBase& mesh,
            EquationSystems* equation_systems,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    file_name = data_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);
    file_name = data_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
} // output_data
