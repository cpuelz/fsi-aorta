// Copyright (c) 2011-2013, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef four_chambered_heart_model_initialization_h
#define four_chambered_heart_model_initialization_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>

// LIBMESH INCLUDES
#include <libmesh/exodusII_io.h>

// other includes
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/PartContext.h>

// SR
// Forward declaration
enum class Active_Strain_Model;
// SR

// ModelInitialization is a static class that provides data and functions
// required to initialize the heart model.
class ModelInitialization
{
public:
    static void transform_meshes(std::vector<MeshBase*>& mesh_vector, 
                                 Pointer<AppInitializer> app_initializer,
                                 const std::vector<std::string>& part_names);

    static void uniformly_refine_meshes(std::vector<MeshBase*>& mesh_vector,
                                        Pointer<AppInitializer> app_initializer,
                                        const std::vector<std::string>& part_names);
    
    static void locally_refine_meshes(std::vector<MeshBase*>& mesh_vector,
                                      Pointer<AppInitializer> app_initializer,
                                      const std::vector<std::string>& part_names);
 
    static void reinit_equation_systems(std::vector<EquationSystems*>& eq_system_vector,
                                        Pointer<AppInitializer> app_initializer,
                                        const std::vector<std::string>& part_names);
    
    static void init_fiber_system(ExodusII_IO& mesh_reader,
                                  EquationSystems& system,
                                  Pointer<Database> input_db,
                                  const std::size_t sys_size);

    // for meshes in exodus format
    static void load_fiber_system(ExodusII_IO& mesh_reader, EquationSystems& system, Pointer<Database> input_db);

    static void load_mesh_info(std::vector<MeshBase*>& mesh_vector,
                               Pointer<AppInitializer> app_initializer,
                               const std::vector<std::string>& part_names);

    static void load_model_parameters(Pointer<AppInitializer> app_initializer,
                                      const std::vector<std::string>& part_names,
                                      std::vector<std::unique_ptr<PartContext> >& context_vector);

    static void create_part_dictionary();

    // SR
    // Active Strain map to set the model from the input file
    typedef std::map<std::string, Active_Strain_Model> ActiveStrainMap;
    static ActiveStrainMap active_strain_model_map;
    // SR
    //static void init_active_strain_map();

private:
    ModelInitialization();
    ModelInitialization(ModelInitialization&);
    ~ModelInitialization();
    ModelInitialization& operator=(ModelInitialization&);

};

#endif // four_chambered_heart_model_initialization_h
