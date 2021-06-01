
/*
 * File:   MeshInfo.cpp
 * Author: cpuelz
 *
 * Created on November 30, 2017, 1:22 PM
 */

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

#include <four_chambered_heart/MeshInfo.h>

// declaration of static variables
std::vector<std::map<std::string, int> > MeshInfo::block_map;
std::vector<std::map<std::string, int> > MeshInfo::sideset_map;
std::vector<std::map<std::string, int> > MeshInfo::nodeset_map;

std::vector<bool> MeshInfo::has_fibers;
std::vector<bool> MeshInfo::turn_on_meters;
std::vector<bool> MeshInfo::turn_on_sources;

std::vector<std::vector<int> > MeshInfo::surface_pressure_IDs;
std::vector<std::vector<int> > MeshInfo::tether_body_force_IDs;
std::vector<std::vector<int> > MeshInfo::circ_model_body_force_IDs;
std::vector<std::vector<int> > MeshInfo::tether_surface_force_IDs;

std::vector<SAMRAI::tbox::Array<std::string> > MeshInfo::tether_body_force_names;
std::vector<SAMRAI::tbox::Array<std::string> > MeshInfo::tether_surface_force_names;
std::vector<SAMRAI::tbox::Array<std::string> > MeshInfo::surface_pressure_names;

bool MeshInfo::circ_model_surface_tethering;
bool MeshInfo::circ_model_body_tethering;
std::map<int, Afterload_Parms> MeshInfo::subdomain_to_afterload_parms;
std::vector<double> MeshInfo::Cartesian_L(NDIM, 0.0);

std::vector<std::vector<int> > MeshInfo::local_refinement_subdomain_IDs;
std::vector<SAMRAI::tbox::Array<std::string> > MeshInfo::local_refinement_subdomain_names;

std::vector<std::vector<std::string> > MeshInfo::subdomain_names;
std::vector<std::vector<libMesh::subdomain_id_type> > MeshInfo::subdomain_IDs;

std::vector<SAMRAI::tbox::Array<std::string> > MeshInfo::afterload_model_nodeset_names;
std::vector<Afterload_Parms> MeshInfo::afterload_parms;

std::vector<int> MeshInfo::pericardial_tethering_IDs;
SAMRAI::tbox::Array<std::string> MeshInfo::pericardial_tethering_names;

int MeshInfo::HEART               = -1000000;
int MeshInfo::AORTA               = -1000000;
int MeshInfo::PULM_ART            = -1000000;
int MeshInfo::AORTA_CAP           = -1000000;
int MeshInfo::PULM_ART_CAP        = -1000000;
int MeshInfo::AORTIC_VALVE        = -1000000;
int MeshInfo::PULM_VALVE          = -1000000;
int MeshInfo::VEIN_CAPS           = -1000000;
int MeshInfo::MITRAL_VALVE        = -1000000;
int MeshInfo::MITRAL_PAP          = -1000000;
int MeshInfo::MITRAL_CHORDS       = -1000000;
int MeshInfo::MITRAL_STRUT_CHORDS = -1000000;
int MeshInfo::TRI_VALVE           = -1000000;
int MeshInfo::TRI_PAP             = -1000000;
int MeshInfo::TRI_CHORDS          = -1000000;
int MeshInfo::PULM_ART_WITH_VALVE = -1000000;
int MeshInfo::AORTA_WITH_VALVE    = -1000000;
int MeshInfo::EVERYTHING_ELSE     = -1000000;

std::map<std::string, int*> MeshInfo::part_dictionary;
