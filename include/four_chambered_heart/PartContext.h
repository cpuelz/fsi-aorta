
/*
 * PartContext.h
 * This class is intended to be used when defining stress functions in
 * order to reduce the number of redundant functions otherwise needed
 * within the code. The class is passed into the stress functions via
 * void* ctx. The part_id corresponds to the index of the part with
 * respect to the parameter vectors. The ModelParms maps carry the
 * pointers to the structs, so the PartContext class is ultimately
 * responsible for the destruction of the structs. Each PartContext
 * carries all of the parameter information for all of the included
 * subdomains. The maps map the subdomain IDs to parameter structs.
 */

#ifndef four_chambered_heart_part_context_h
#define four_chambered_heart_part_context_h

#include <four_chambered_heart/ModelParameters.h>

class PartContext
{
public:
    explicit PartContext(int id) : part_id(id)
    {
    }

    // part index
    int part_id;

    // Constitutive model ID map
    std::unordered_map<libMesh::subdomain_id_type, int> constitutive_model_ID;

    // Constitutive model ID map
    std::unordered_map<libMesh::subdomain_id_type, std::string> constitutive_model_name;

    // dilational stress ID map
    std::unordered_map<libMesh::subdomain_id_type, int> volumetric_energy_ID;

    // General model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<General_ModelParms> > general_parms;

    // Guccione model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<Guccione_ModelParms> > guccione_parms;

    // HGO model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<HGO_ModelParms> > hgo_parms;

    // HO model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<HO_ModelParms> > ho_parms;

    // Neohookean model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<Neohookean_ModelParms> > neohookean_parms;

    // Exponential neohookean model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<Exp_Neohookean_ModelParms> > exp_neohookean_parms;

    // Biaxial fung type model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<Biaxial_Fung_Type_ModelParms> >
        biaxial_fung_type_parms;

    // Augustin model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<Augustin_ModelParms> > augustin_parms;

    // nonlinear spring model parameter map
    std::unordered_map<libMesh::subdomain_id_type, std::unique_ptr<NonlinearSpring_ModelParms> > nonlinear_spring_parms;

    // store info about invariants in each subdomain of each part
    std::unordered_map<libMesh::subdomain_id_type, double > min_J;
    std::unordered_map<libMesh::subdomain_id_type, double > max_J;
    std::unordered_map<libMesh::subdomain_id_type, double > min_I1;
    std::unordered_map<libMesh::subdomain_id_type, double > max_I1;
    std::unordered_map<libMesh::subdomain_id_type, double > min_I4_one;
    std::unordered_map<libMesh::subdomain_id_type, double > max_I4_one;
    std::unordered_map<libMesh::subdomain_id_type, double > min_I4_two;
    std::unordered_map<libMesh::subdomain_id_type, double > max_I4_two;
    std::unordered_map<libMesh::subdomain_id_type, double > min_I4_three;
    std::unordered_map<libMesh::subdomain_id_type, double > max_I4_three;

    // store active tension data for each subdomain of each part  
    std::unordered_map<libMesh::subdomain_id_type, double > active_tension;

private:
};

#endif
