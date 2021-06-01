/*
 * File:   ModelParameters.h
 * Author: cpuelz
 *
 * Created on April 23, 2018, 2:28 PM
 */

#ifndef four_chambered_heart_model_parameters_h
#define four_chambered_heart_model_parameters_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

#include <libmesh/point.h>
#include <libmesh/elem.h>
#include <string>

struct Original_ActiveTensionParms
{
    double t_active;
    double t_relax;
    double t_ramp;
    double t_period;
    double t_delay;
};

struct Augustin_ActiveTensionParms
{
    double t_ramp;
    double t_period;
    double t_delay;
    double t_t;
    double tau_c;
    double tau_r;
};

struct Data_Activation_FunctionParms
{
    double t_active;
    double t_ramp;
    double t_period;
    double t_delay;
};

struct Pressure_Activation_FunctionParms
{
    double t_peak;
    double t_drop;
    double t_plateau;
    double t_ramp;
    double t_period;
    double t_delay;
};

// read in the input file as in the comments
enum class Active_Strain_Model { TRANSVERSELY_ISOTROPIC,   // transverse
                                 ORTHOTROPIC,              // orthotropic
                                 TRANSMURALLY_ORTHOTROPIC};// transmural

struct General_ModelParms
{
    // for volumetric energy
    double beta_s;

    // for prestrain (should probably switch this to the HO model parms)
    double t_ramp;

    // parms for active tension
    double Tension;
    Original_ActiveTensionParms original_active_tension_parms;
    Augustin_ActiveTensionParms augustin_active_tension_parms;
    Data_Activation_FunctionParms data_activation_function_parms;
    Pressure_Activation_FunctionParms pressure_activation_function_parms;
    unsigned int first_active_tension_fiber_index;
    unsigned int second_active_tension_fiber_index;
    unsigned int third_active_tension_fiber_index;
    double (*active_tension_function)(double time,
                                      Elem* const elem,
                                      void* ctx);

    // ACTIVE STRESS PARAMETERS
    bool enable_active_stress;
    bool scaling_with_I4;
    double nu1;
    double nu2;
    double nu3;

    // ACTIVE STRAIN PARAMETERS
    bool enable_active_strain;
    // scale activation function by this number
    double max_fiber_shortening;
    // Link between the microscale and macroscale for the orthotropic models
    double active_strain_kappa;
    // choose between the models
    Active_Strain_Model active_strain_model;
};


struct HGO_ModelParms
{
    double c;
    double c0;
    double c2;
    unsigned int radial_fiber_index;
    unsigned int circ_fiber_index;
    double radial_prestrain_coeff;
    double circ_prestrain_coeff;
};

struct HO_ModelParms
{
    double mu_e;
    double a, b;
    double af, bf;
    double as, bs;
    double afs, bfs;
    double kf, ks; // fiber dispersion parameters
    bool single_fiber_family;
    bool use_scaled_invariants;
    bool turn_off_fibers_in_compression;
    bool DEV_projection;
    bool enable_prestrain;
    unsigned int radial_fiber_index;
    unsigned int circ_fiber_index;
    double radial_prestrain_coeff;
    double circ_prestrain_coeff;
};

struct Guccione_ModelParms
{
    double bf;
    double bt;
    double bfs;
    double C;
    double circ_prestrain_coeff;
    unsigned int circ_fiber_index;
    bool DEV_projection;
};

struct Neohookean_ModelParms
{
    double mu_s;
};

struct Exp_Neohookean_ModelParms
{
    double a;
    double b;
};

struct Biaxial_Fung_Type_ModelParms
{
    double C;
    double A1;
    double A2;
    double A3;
    double A4;
    double A5;
    double A6;
    unsigned int radial_fiber_index;
    unsigned int circ_fiber_index;
    bool DEV_projection;
};

struct Augustin_ModelParms
{
    double kappa;
    double a;
    double b;
    double af;
    double bf;
    int circ_fiber_index;
    bool DEV_projection;
};

struct NonlinearSpring_ModelParms
{
    double af;
    double bf;
    int fiber_index;
    double prestrain_coeff;
    bool enable_prestrain;
    bool use_scaled_invariants;
    bool DEV_projection;
};

struct Afterload_Parms
{
    std::string nodeset_name;
    int nodeset_ID;
    libMesh::Point centroid;
    double radius;
    int axis;
    int side;
    int part;
};

#endif // four_chambered_heart_model_parameters_h
