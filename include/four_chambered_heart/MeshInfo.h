
/*
 * File:   MeshInfo.h
 * Author: cpuelz
 *
 * Created on November 30, 2017, 1:22 PM
 */

#ifndef four_chambered_heart_mesh_info_h
#define four_chambered_heart_mesh_info_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

#include <map>
#include <stddef.h>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>
#include <libmesh/point.h>
#include <tbox/Array.h>

#include <four_chambered_heart/ModelParameters.h>

class MeshInfo
{
public:
    
    // maps to store info
    static std::vector<std::map<std::string, int> > block_map;
    static std::vector<std::map<std::string, int> > sideset_map;
    static std::vector<std::map<std::string, int> > nodeset_map;

    // vector to store boolean for determining if a given mesh
    // has fiber systems stored/
    static std::vector<bool> has_fibers;
    static std::vector<bool> turn_on_meters;
    static std::vector<bool> turn_on_sources;

    // vectors to store subdomain and sideset IDs for prescribing boundary conditions,
    // active contraction, and tethering forces
    static std::vector<std::vector<int> > surface_pressure_IDs;
    static std::vector<std::vector<int> > tether_body_force_IDs;
    static std::vector<std::vector<int> > circ_model_body_force_IDs;
    static std::vector<std::vector<int> > tether_surface_force_IDs;
    
    // whether to do surface and/or body tethering with the circulation models
    static bool circ_model_surface_tethering;
    static bool circ_model_body_tethering;
    static std::map<int, Afterload_Parms> subdomain_to_afterload_parms;
    static std::vector<double> Cartesian_L;
    
    // arrays to store names of different mesh sidesets and blocksets for body forces and boundary conditions
    static std::vector<SAMRAI::tbox::Array<std::string> > tether_body_force_names;
    static std::vector<SAMRAI::tbox::Array<std::string> > tether_surface_force_names;
    static std::vector<SAMRAI::tbox::Array<std::string> > surface_pressure_names;
    
    // subdomain names and IDs for local mesh refinement
    static std::vector<std::vector<int> > local_refinement_subdomain_IDs;
    static std::vector<SAMRAI::tbox::Array<std::string> > local_refinement_subdomain_names;

    // for subdomain names
    static std::vector<std::vector<std::string> > subdomain_names;

    // for subdomain IDs
    static std::vector<std::vector<libMesh::subdomain_id_type> > subdomain_IDs;

    // for afterload model names
    static std::vector<SAMRAI::tbox::Array<std::string> > afterload_model_nodeset_names;

    // a vector of afterload model parm structs
    static std::vector<Afterload_Parms> afterload_parms;

    // these are used only for the HeartPart
    static std::vector<int> pericardial_tethering_IDs;
    static SAMRAI::tbox::Array<std::string> pericardial_tethering_names;

    // integer IDs for different parts
    static int HEART;
    static int AORTA;
    static int PULM_ART;
    static int AORTA_CAP;
    static int PULM_ART_CAP;
    static int AORTIC_VALVE;
    static int PULM_VALVE;
    static int VEIN_CAPS;
    static int MITRAL_VALVE;
    static int MITRAL_PAP;
    static int MITRAL_CHORDS;
    static int MITRAL_STRUT_CHORDS;
    static int TRI_VALVE;
    static int TRI_PAP;
    static int TRI_CHORDS;
    static int PULM_ART_WITH_VALVE;
    static int AORTA_WITH_VALVE;
    static int EVERYTHING_ELSE;

    static std::map<std::string, int*> part_dictionary;

    MeshInfo();
    MeshInfo(const MeshInfo& orig);
    virtual ~MeshInfo();

private:
};

#endif // four_chambered_heart_mesh_info_h
