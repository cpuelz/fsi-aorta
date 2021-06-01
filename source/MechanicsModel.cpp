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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>

// APPLICATION INCLUDES
#include <four_chambered_heart/BoundaryConditions.h>
#include <four_chambered_heart/MechanicsModel.h>
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/ModelInitialization.h>
#include <four_chambered_heart/ModelParameters.h>
#include <four_chambered_heart/PartContext.h>

#include <cmath>

// CLASS IMPLEMENTATION
std::string MechanicsModel::HO_NAME = "ho_model";
std::string MechanicsModel::HGO_NAME = "hgo_model";
std::string MechanicsModel::GUCCIONE_NAME = "guccione_model";
std::string MechanicsModel::NEOHOOKEAN_NAME = "neohookean_model";
std::string MechanicsModel::EXP_NEOHOOKEAN_NAME = "exp_neohookean_model";
std::string MechanicsModel::BIAXIAL_FUNG_TYPE_NAME = "biaxial_fung_type_model";
std::string MechanicsModel::AUGUSTIN_NAME = "augustin_model";
std::string MechanicsModel::NONLINEAR_SPRING_NAME = "nonlinear_spring_model";

int MechanicsModel::HO_ID = 0;
int MechanicsModel::HGO_ID = 1;
int MechanicsModel::GUCCIONE_ID = 2;
int MechanicsModel::NEOHOOKEAN_ID = 3;
int MechanicsModel::EXP_NEOHOOKEAN_ID = 4;
int MechanicsModel::BIAXIAL_FUNG_TYPE_ID = 5;
int MechanicsModel::AUGUSTIN_ID = 6;
int MechanicsModel::NONLINEAR_SPRING_ID = 7;

std::string MechanicsModel::ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME = "original_active_tension";
std::string MechanicsModel::AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME = "augustin_active_tension";
std::string MechanicsModel::DATA_ACTIVE_TENSION_FUNCTION_NAME = "data_active_tension";
std::string MechanicsModel::PRESSURE_ACTIVE_TENSION_FUNCTION_NAME = "pressure_active_tension";

std::string MechanicsModel::LOG_VOLUMETRIC_ENERGY_NAME = "log_volumetric_energy";
std::string MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_NAME = "quadratic_volumetric_energy";

int MechanicsModel::LOG_VOLUMETRIC_ENERGY_ID = 0;
int MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_ID = 1;

bool MechanicsModel::inverted_element = false;

void
MechanicsModel::get_PK1_dev_stress_function_systems(vector<unsigned int>& systems)
{
    systems.resize(0);
    return;
}

double
MechanicsModel::original_active_tension_function(const double time,
                                                 Elem* const elem,
                                                 void* ctx)
{
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const double t_active = part_context->general_parms.at(subdomain_ID)->original_active_tension_parms.t_active;
    const double t_period = part_context->general_parms.at(subdomain_ID)->original_active_tension_parms.t_period;
    const double t_relax = part_context->general_parms.at(subdomain_ID)->original_active_tension_parms.t_relax;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->original_active_tension_parms.t_ramp;
    const double t_delay = part_context->general_parms.at(subdomain_ID)->original_active_tension_parms.t_delay;
    const double rem = fmod(time - t_ramp, t_period);

    if (time < t_ramp + t_delay)
    {
        return 0.0;
    }
    else if (rem > t_delay && rem < t_active + t_relax + t_delay)
    {
        if (rem < t_active + t_delay)
        {
            return (rem - t_delay) / t_active;
        }
        else
        {
            return (t_active + t_relax + t_delay - rem) / (t_relax);
        }
    }
    else
    {
        return 0.0;
    }
}

double
MechanicsModel::augustin_active_tension_function(const double time,
                                                 Elem* const elem,
                                                 void* ctx)
{
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const double t_t = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.t_t;
    const double tau_c = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.tau_c;
    const double tau_r = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.tau_r;
    const double t_period = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.t_period;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.t_ramp;
    const double t_delay = part_context->general_parms.at(subdomain_ID)->augustin_active_tension_parms.t_delay;
    const double rem = fmod(time - t_ramp, t_period);

    if (time < t_ramp + t_delay)
    {
        return 0.0;
    }
    else if (rem > t_delay && rem < t_t + t_delay)
    {
        //return Tension * tanh(rem/tau_c) * tanh(rem/tau_c) * tanh((t_t - rem)/tau_r) * tanh((t_t - rem)/tau_r);
        double tanh1 = std::tanh(rem/tau_c);
        double tanh2 = std::tanh((t_t - rem)/tau_r);
        return tanh1 * tanh1 * tanh2 * tanh2;
    }
    else
    {
        return 0.0;
    }
}



// Activation function taken from experimental data
// The activation function is scaled over a time unit
// To scale it: divide the time by the desired period
static std::array<double, 24> activation_times = {
        0,
        0.0330,
        0.0413,
        0.0606,
        0.0771,
        0.1019,
        0.1211,
        0.1569,
        0.1817,
        0.2065,
        0.2423,
        0.2643,
        0.2781,
        0.2863,
        0.3001,
        0.3221,
        0.3579,
        0.4350,
        0.5176,
        0.6195,
        0.7213,
        0.8177,
        0.9333,
        1.0000
};

// magnitude of the activation
static std::array<double, 24> activation = {
        0,
        0.0392,
        0.1047,
        0.2233,
        0.4024,
        0.5739,
        0.7328,
        0.8740,
        0.9497,
        0.9799,
        1.0000,
        0.9747,
        0.9469,
        0.8662,
        0.7526,
        0.6365,
        0.4699,
        0.2957,
        0.1921,
        0.1086,
        0.0453,
        0.0224,
        0.0086,
        0
};


double
MechanicsModel::data_activation_function ( const double time,
                                            Elem* const elem,
                                            void* ctx)
{
    // Get parameters
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const double t_active = part_context->general_parms.at(subdomain_ID)->data_activation_function_parms.t_active;
    const double t_period = part_context->general_parms.at(subdomain_ID)->data_activation_function_parms.t_period;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->data_activation_function_parms.t_ramp;
    const double t_delay = part_context->general_parms.at(subdomain_ID)->data_activation_function_parms.t_delay;

    // if we are in the ramp phase do not contract
    if (time < t_ramp + t_delay) return 0.0;
    // get time with respect to period
    double t = fmod(time-t_ramp-t_delay, t_period);
    // if we are in the relaxation phase do not contract
    if (t > t_active )  return 0.0;

    // at what percentage of the contraction time are we?
    double scaled_t = t / t_active;

    // get iterator to the first data point larger than the current scaled time
    auto it_t1=std::lower_bound(activation_times.begin(), activation_times.end(), scaled_t); //
    // get iterator to the previous data point
    auto it_t0=it_t1-1;

    // find the index in the array corresponding to the iterators
    int i1 = std::distance(activation_times.begin(), it_t1);
    int i0 = i1-1;

    // Get the time interval in which the scaled_t falls
    double t0 = *it_t0;
    double t1 = *it_t1;

    // get the (normalized) force interval
    double f0 = activation[i0];
    double f1 = activation[i1];

    // interpolate linearly
    double f =  (scaled_t - t0) * (f1 - f0) / (t1 - t0) + f0;

    return f;
}


double
MechanicsModel::pressure_activation_function ( const double time,
                                            Elem* const elem,
                                            void* ctx)
{
    // Get parameters
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const double t_period = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_period;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_ramp;
    const double t_delay = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_delay;
    const double t_peak = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_peak;
    const double t_plateau = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_plateau;
    const double t_drop = part_context->general_parms.at(subdomain_ID)->pressure_activation_function_parms.t_drop;

    // if we are in the ramp phase do not contract
    if (time < t_ramp + t_delay) return 0.0;
    // get time with respect to period
    double t = fmod(time-t_ramp-t_delay, t_period);
    // if we are in the relaxation phase do not contract
    if (t >= t_peak + t_drop + t_plateau )  return 0.0;

    double f = 0.0;
    if( t <= t_peak ) f = ( 0.5 - 0.5 * std::cos( M_PI * t / t_peak) );
    else if ( t >=  t_peak + t_plateau ) f = ( 0.5 + 0.5 * std::cos( M_PI * t / t_drop - M_PI / t_drop * ( t_plateau + t_peak ) ) );
    else f = 1.0;

    return f;
}

void
MechanicsModel::PK1_dev_stress_function(TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& x, // current location
                                        const libMesh::Point& X, // reference location
                                        Elem* const elem,
                                        const std::vector<const std::vector<double>*>& var_data,
                                        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                        double data_time,
                                        void* ctx)
{
    PP.zero();
    TensorValue<double> PP_active; PP_active.zero();
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const int part_ID = part_context->part_id;
    const int constitutive_model_ID = part_context->constitutive_model_ID.at(subdomain_ID);
    const bool enable_active_stress = part_context->general_parms.at(subdomain_ID)->enable_active_stress;
    const bool enable_active_strain = part_context->general_parms.at(subdomain_ID)->enable_active_strain;
    const unsigned int first_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->first_active_tension_fiber_index;
    const unsigned int second_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->second_active_tension_fiber_index;
    const unsigned int third_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->third_active_tension_fiber_index;
    const TensorValue<double> CC = FF.transpose() * FF;

    // Common setup for active stress and active strain
    // Get the fiber directions
    libMesh::VectorValue<double> f0;
    libMesh::VectorValue<double> s0;
    libMesh::VectorValue<double> n0;

    libMesh::VectorValue<double> f;
    libMesh::VectorValue<double> s;
    libMesh::VectorValue<double> n;

    // Set up models
    double normalized_activation = 0.0;
    if(enable_active_strain || enable_active_stress)
    {
        TBOX_ASSERT(0 < var_data.size());
        TBOX_ASSERT(9 <= var_data[0]->size());

        // active stress in "main" fiber direction
        f0(0) = (*var_data[0])[0 + NDIM * first_active_tension_fiber_index];
        f0(1) = (*var_data[0])[1 + NDIM * first_active_tension_fiber_index];
        f0(2) = (*var_data[0])[2 + NDIM * first_active_tension_fiber_index];
        f0 = f0.unit();
        f = FF * f0;

        s0(0) = (*var_data[0])[0 + NDIM * second_active_tension_fiber_index];
        s0(1) = (*var_data[0])[1 + NDIM * second_active_tension_fiber_index];
        s0(2) = (*var_data[0])[2 + NDIM * second_active_tension_fiber_index];
        s0 = s0.unit();
        s = FF * s0;

        n0(0) = (*var_data[0])[0 + NDIM * third_active_tension_fiber_index];
        n0(1) = (*var_data[0])[1 + NDIM * third_active_tension_fiber_index];
        n0(2) = (*var_data[0])[2 + NDIM * third_active_tension_fiber_index];
        n0 = n0.unit();
        n = FF * n0;

        // evaluate activation function
        normalized_activation = part_context->general_parms.at(subdomain_ID)
                ->active_tension_function(data_time, elem, ctx);

        // logging I4n
        part_context->min_I4_three.at(subdomain_ID) = std::min(part_context->min_I4_three.at(subdomain_ID), n * n);
        part_context->max_I4_three.at(subdomain_ID) = std::max(part_context->max_I4_three.at(subdomain_ID), n * n);
    }

    // stuff for logging min/max of the invariants
    part_context->min_J.at(subdomain_ID) = std::min(part_context->min_J.at(subdomain_ID), FF.det());
    part_context->max_J.at(subdomain_ID) = std::max(part_context->max_J.at(subdomain_ID), FF.det());
    part_context->min_I1.at(subdomain_ID) = std::min(part_context->min_I1.at(subdomain_ID), CC.tr());
    part_context->max_I1.at(subdomain_ID) = std::max(part_context->max_I1.at(subdomain_ID), CC.tr());

    // Active Strain setup:
    // FF = FFE * FFA
    // Evaluate PPE = PPE(FFE) using FFE = FF * FFAinv
    // and then pull back PP = PPE * FAinv
    libMesh::TensorValue<double> FFA(1,0,0,
                                     0,1,0,
                                     0,0,1);
    libMesh::TensorValue<double> FFAinv;
    //  setup FFA
    if (enable_active_strain)
    {
        // scale activation function by this number
        const double max_fiber_shortening   = part_context->general_parms.at(subdomain_ID)->max_fiber_shortening;
        // Link between the microscale and macroscale for the orthotropic models
        const double kappa   = part_context->general_parms.at(subdomain_ID)->active_strain_kappa;
        // choose between the models
        const Active_Strain_Model active_strain_model  = part_context->general_parms.at(subdomain_ID)->active_strain_model;
        // get parameter for the transmural model (Currently disabled)
        // TO DO:
        // Enable this
        const double lambda = 0.0;
        // const double lambda = (*var_data[0])[9];
        // setup FFA
        double gamma_f = -max_fiber_shortening * normalized_activation;
        double gamma_s = 0.0;
        double gamma_n = 0.0;
        switch(active_strain_model)
        {
            case Active_Strain_Model::ORTHOTROPIC:
            {
                gamma_n = kappa * gamma_f;
                gamma_s = 1.0 / (1.0 + gamma_f) / (1.0 + gamma_n) - 1.0;
                break;
            }
            case Active_Strain_Model::TRANSMURALLY_ORTHOTROPIC:
            {
                // lambda = 0 on the endocardium
                // lambda = 1 on the epicardium
                // there is a typo in the paper ... Luca Luca Luca ...
                gamma_n = (1.0 - lambda) * kappa * gamma_f
                        + lambda * ( 1.0 / std::sqrt( 1 + gamma_f ) - 1.0 );
                gamma_s = 1.0 / (1.0 + gamma_f) / (1.0 + gamma_n) - 1.0;
                break;
            }
            case Active_Strain_Model::TRANSVERSELY_ISOTROPIC:
            default:
            {
                gamma_s = 1.0 / std::sqrt(1.0 + gamma_f) - 1.0;
                gamma_n = gamma_s;
                break;
            }
        }
        FFA += gamma_f * outer_product(f0, f0);
        FFA += gamma_s * outer_product(s0, s0);
        FFA += gamma_n * outer_product(n0, n0);
        FFAinv = FFA.inverse();
    }

    // Hack to reduce the amount of copy and paste of functions in the code:
    // Elastic Stress
    libMesh::TensorValue<double> PPE(PP);
    // Elastic Defomration Gradient tensor
    libMesh::TensorValue<double> FFE(FF);
    if (enable_active_strain) FFE = FFE * FFAinv;

    if (constitutive_model_ID == MechanicsModel::HO_ID)
    {
        PK1_dev_stress_function_general_HO(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::HGO_ID)
    {
        // need to implement
    }
    if (constitutive_model_ID == MechanicsModel::GUCCIONE_ID)
    {
        PK1_dev_stress_function_general_Guccione(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::NEOHOOKEAN_ID)
    {
        PK1_dev_stress_function_general_neohookean(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::EXP_NEOHOOKEAN_ID)
    {
        PK1_dev_stress_function_general_exp_neohookean(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::BIAXIAL_FUNG_TYPE_ID)
    {
        PK1_dev_stress_function_general_biaxial_fung_type(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::AUGUSTIN_ID)
    {
        PK1_dev_stress_function_general_augustin(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::NONLINEAR_SPRING_ID)
    {
        PK1_dev_stress_function_general_nonlinear_spring(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }

    if (enable_active_strain) PP += PPE * FFAinv;
    else PP = PPE;
    // finished hack

    // Add active stress
    if (enable_active_stress)
    {
        // this is Tmax * normalized_function_of_time
        const double Tension = normalized_activation * part_context->general_parms.at(subdomain_ID)->Tension;
        part_context->active_tension.at(subdomain_ID) = std::max(part_context->active_tension.at(subdomain_ID), Tension);
        // percentage contribution in the various directions
        // should be in [0, 1]
        const double nu1 = part_context->general_parms.at(subdomain_ID)->nu1;
        const double nu2 = part_context->general_parms.at(subdomain_ID)->nu2;
        const double nu3 = part_context->general_parms.at(subdomain_ID)->nu3;
        const bool scaling_with_I4 = part_context->general_parms.at(subdomain_ID)->scaling_with_I4;

        const double I4f = scaling_with_I4 ? f.norm() : 1.0;
        PP_active += Tension * ( nu1 / sqrt(I4f) ) * outer_product(f, f0);

        // contribution of active stress in other fiber directions
        if (!MathUtilities<double>::equalEps(nu2, 0.0))
        {
            const double I4s = scaling_with_I4 ? s.norm() : 1.0;
            PP_active += Tension * ( nu2 / sqrt(I4s) ) * outer_product(s, s0);
        }

        if (!MathUtilities<double>::equalEps(nu3, 0.0))
        {
            const double I4n = scaling_with_I4 ? n.norm() : 1.0;
            PP_active += Tension * ( nu3 / sqrt(I4n) ) * outer_product(n, n0);
        }

        // compute deviatoric projection of the active stress
        compute_deviatoric_projection(PP_active, FF);
        PP += PP_active;
    }
    if (FF.det() <= 0.01)
    {
        inverted_element = true;
        PP.zero();
        TBOX_WARNING("an element is close to inverting, or "
                     "may have inverted, in subdomain "
                     << subdomain_ID << " of part "
                     << part_ID << " !!!!!\n");
    }
    return;
}

void
MechanicsModel::PK1_dil_stress_function(TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& /*x*/, // current location
                                        const libMesh::Point& /*X*/, // reference location
                                        Elem* const elem,
                                        const std::vector<const std::vector<double>*>& /*var_data*/,
                                        const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                                        double /*data_time*/,
                                        void* ctx)
{
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const double J = FF.det();
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();
    const double beta_s = part_context->general_parms.at(subdomain_ID)->beta_s;
    const int volumetric_energy_ID = part_context->volumetric_energy_ID.at(subdomain_ID);

    PP.zero();
    if (volumetric_energy_ID == MechanicsModel::LOG_VOLUMETRIC_ENERGY_ID)
    {
        if (!MathUtilities<double>::equalEps(beta_s, 0.0))
        {
            // W(J) = beta_s*(J * log(J) - J + 1)
            PP += beta_s * J * log(J) * FF_inv_trans;
        }
    }
    if (volumetric_energy_ID == MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_ID)
    {
       if (!MathUtilities<double>::equalEps(beta_s, 0.0))
        {
            // W(J) = 0.5 * beta_s * (J - 1)^2
            PP += beta_s * J * (J - 1.0) * FF_inv_trans;
        }
    }

    // check to see if J is small.
    if (J <= 0.01)
    {
        inverted_element = true;
        PP.zero();
    }
    return;
}

void
MechanicsModel::PK1_combined_stress_function(TensorValue<double>& PP,
                                             const TensorValue<double>& FF,
                                             const libMesh::Point& x, // current location
                                             const libMesh::Point& X, // reference location
                                             Elem* const elem,
                                             const std::vector<const std::vector<double>*>& var_data,
                                             const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                             double data_time,
                                             void* ctx)
{
    PP.zero();
    TensorValue<double> PP_active; PP_active.zero();
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();
    const int part_ID = part_context->part_id;
    const int constitutive_model_ID = part_context->constitutive_model_ID.at(subdomain_ID);
    const bool enable_active_stress = part_context->general_parms.at(subdomain_ID)->enable_active_stress;
    const bool enable_active_strain = part_context->general_parms.at(subdomain_ID)->enable_active_strain;
    const unsigned int first_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->first_active_tension_fiber_index;
    const unsigned int second_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->second_active_tension_fiber_index;
    const unsigned int third_active_tension_fiber_index =
        part_context->general_parms.at(subdomain_ID)->third_active_tension_fiber_index;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const double J = FF.det();
    const TensorValue<double> CC = FF.transpose() * FF;
    const double beta_s = part_context->general_parms.at(subdomain_ID)->beta_s;
    const int volumetric_energy_ID = part_context->volumetric_energy_ID.at(subdomain_ID);

    // Common setup for active stress and active strain
    // Get the fiber directions
    libMesh::VectorValue<double> f0;
    libMesh::VectorValue<double> s0;
    libMesh::VectorValue<double> n0;

    libMesh::VectorValue<double> f;
    libMesh::VectorValue<double> s;
    libMesh::VectorValue<double> n;

    // Set up models
    double normalized_activation = 0.0;
    if(enable_active_strain || enable_active_stress)
    {
        TBOX_ASSERT(0 < var_data.size());
        TBOX_ASSERT(9 <= var_data[0]->size());

        // active stress in "main" fiber direction
        f0(0) = (*var_data[0])[0 + NDIM * first_active_tension_fiber_index];
        f0(1) = (*var_data[0])[1 + NDIM * first_active_tension_fiber_index];
        f0(2) = (*var_data[0])[2 + NDIM * first_active_tension_fiber_index];
        f0 = f0.unit();
        f = FF * f0;

        s0(0) = (*var_data[0])[0 + NDIM * second_active_tension_fiber_index];
        s0(1) = (*var_data[0])[1 + NDIM * second_active_tension_fiber_index];
        s0(2) = (*var_data[0])[2 + NDIM * second_active_tension_fiber_index];
        s0 = s0.unit();
        s = FF * s0;

        n0(0) = (*var_data[0])[0 + NDIM * third_active_tension_fiber_index];
        n0(1) = (*var_data[0])[1 + NDIM * third_active_tension_fiber_index];
        n0(2) = (*var_data[0])[2 + NDIM * third_active_tension_fiber_index];
        n0 = n0.unit();
        n = FF * n0;

        // evaluate activation function
        normalized_activation = part_context->general_parms.at(subdomain_ID)
                ->active_tension_function(data_time, elem, ctx);

        // logging I4n
        part_context->min_I4_three.at(subdomain_ID) = std::min(part_context->min_I4_three.at(subdomain_ID), n * n);
        part_context->max_I4_three.at(subdomain_ID) = std::max(part_context->max_I4_three.at(subdomain_ID), n * n);
    }

    // stuff for logging min/max of the invariants
    part_context->min_J.at(subdomain_ID) = std::min(part_context->min_J.at(subdomain_ID), FF.det());
    part_context->max_J.at(subdomain_ID) = std::max(part_context->max_J.at(subdomain_ID), FF.det());
    part_context->min_I1.at(subdomain_ID) = std::min(part_context->min_I1.at(subdomain_ID), CC.tr());
    part_context->max_I1.at(subdomain_ID) = std::max(part_context->max_I1.at(subdomain_ID), CC.tr());

    // Active Strain setup:
    // FF = FFE * FFA
    // Evaluate PPE = PPE(FFE) using FFE = FF * FFAinv
    // and then pull back PP = PPE * FAinv
    libMesh::TensorValue<double> FFA(1,0,0,
                                     0,1,0,
                                     0,0,1);
    libMesh::TensorValue<double> FFAinv;
    //  setup FFA
    if (enable_active_strain)
    {
        // scale activation function by this number
        const double max_fiber_shortening   = part_context->general_parms.at(subdomain_ID)->max_fiber_shortening;
        // Link between the microscale and macroscale for the orthotropic models
        const double kappa   = part_context->general_parms.at(subdomain_ID)->active_strain_kappa;
        // choose between the models
        const Active_Strain_Model active_strain_model  = part_context->general_parms.at(subdomain_ID)->active_strain_model;
        // get parameter for the transmural model (Currently disabled)
        // TO DO:
        // Enable this
        const double lambda = 0.0;
        // const double lambda = (*var_data[0])[9];
        // setup FFA
        double gamma_f = -max_fiber_shortening * normalized_activation;
        double gamma_s = 0.0;
        double gamma_n = 0.0;
        switch(active_strain_model)
        {
            case Active_Strain_Model::ORTHOTROPIC:
            {
                gamma_n = kappa * gamma_f;
                gamma_s = 1.0 / (1.0 + gamma_f) / (1.0 + gamma_n) - 1.0;
                break;
            }
            case Active_Strain_Model::TRANSMURALLY_ORTHOTROPIC:
            {
                // lambda = 0 on the endocardium
                // lambda = 1 on the epicardium
                // there is a typo in the paper ... Luca Luca Luca ...
                gamma_n = (1.0 - lambda) * kappa * gamma_f
                        + lambda * ( 1.0 / std::sqrt( 1 + gamma_f ) - 1.0 );
                gamma_s = 1.0 / (1.0 + gamma_f) / (1.0 + gamma_n) - 1.0;
                break;
            }
            case Active_Strain_Model::TRANSVERSELY_ISOTROPIC:
            default:
            {
                gamma_s = 1.0 / std::sqrt(1.0 + gamma_f) - 1.0;
                gamma_n = gamma_s;
                break;
            }
        }
        FFA += gamma_f * outer_product(f0, f0);
        FFA += gamma_s * outer_product(s0, s0);
        FFA += gamma_n * outer_product(n0, n0);
        FFAinv = FFA.inverse();
    }

    // Hack to reduce the amount of copy and paste of functions in the code:
    // Elastic Stress
    libMesh::TensorValue<double> PPE(PP);
    // Elastic Defomration Gradient tensor
    libMesh::TensorValue<double> FFE(FF);
    if (enable_active_strain) FFE = FFE * FFAinv;

    if (constitutive_model_ID == MechanicsModel::HO_ID)
    {
        PK1_dev_stress_function_general_HO(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::HGO_ID)
    {
        // need to implement
    }
    if (constitutive_model_ID == MechanicsModel::GUCCIONE_ID)
    {
        PK1_dev_stress_function_general_Guccione(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::NEOHOOKEAN_ID)
    {
        PK1_dev_stress_function_general_neohookean(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::EXP_NEOHOOKEAN_ID)
    {
        PK1_dev_stress_function_general_exp_neohookean(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::BIAXIAL_FUNG_TYPE_ID)
    {
        PK1_dev_stress_function_general_biaxial_fung_type(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::AUGUSTIN_ID)
    {
        PK1_dev_stress_function_general_augustin(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }
    if (constitutive_model_ID == MechanicsModel::NONLINEAR_SPRING_ID)
    {
        PK1_dev_stress_function_general_nonlinear_spring(PPE, FFE, x, X, elem, var_data, grad_var_data, data_time, ctx);
    }

    if (enable_active_strain) PP += PPE * FFAinv;
    else PP = PPE;
    // finished hack

    // Add active stress
    if (enable_active_stress)
    {
        // this is Tmax * normalized_function_of_time
        const double Tension = normalized_activation * part_context->general_parms.at(subdomain_ID)->Tension;
        part_context->active_tension.at(subdomain_ID) = std::max(part_context->active_tension.at(subdomain_ID), Tension);
        // percentage contribution in the various directions
        // should be in [0, 1]
        const double nu1 = part_context->general_parms.at(subdomain_ID)->nu1;
        const double nu2 = part_context->general_parms.at(subdomain_ID)->nu2;
        const double nu3 = part_context->general_parms.at(subdomain_ID)->nu3;
        const bool scaling_with_I4 = part_context->general_parms.at(subdomain_ID)->scaling_with_I4;

        const double I4f = scaling_with_I4 ? f.norm() : 1.0;
        PP_active += Tension * ( nu1 / sqrt(I4f) ) * outer_product(f, f0);

        // contribution of active stress in other fiber directions
        if (!MathUtilities<double>::equalEps(nu2, 0.0))
        {
            const double I4s = scaling_with_I4 ? s.norm() : 1.0;
            PP_active += Tension * ( nu2 / sqrt(I4s) ) * outer_product(s, s0);
        }

        if (!MathUtilities<double>::equalEps(nu3, 0.0))
        {
            const double I4n = scaling_with_I4 ? n.norm() : 1.0;
            PP_active += Tension * ( nu3 / sqrt(I4n) ) * outer_product(n, n0);
        }

        // compute deviatoric projection of the active stress
        compute_deviatoric_projection(PP_active, FF);
        PP += PP_active;
    }
    if (FF.det() <= 0.01)
    {
        inverted_element = true;
        PP.zero();
        TBOX_WARNING("an element is close to inverting, or "
                     "may have inverted, in subdomain "
                     << subdomain_ID << " of part "
                     << part_ID << " !!!!!\n");
    }

    // add dilational stress
    if (volumetric_energy_ID == MechanicsModel::LOG_VOLUMETRIC_ENERGY_ID)
    {
        if (!MathUtilities<double>::equalEps(beta_s, 0.0))
        {
            // W(J) = beta_s*(J * log(J) - J + 1)
            PP += beta_s * J * log(J) * FF_inv_trans;
        }
    }
    if (volumetric_energy_ID == MechanicsModel::QUADRATIC_VOLUMETRIC_ENERGY_ID)
    {
       if (!MathUtilities<double>::equalEps(beta_s, 0.0))
        {
            // W(J) = 0.5 * beta_s * (J - 1)^2
            PP += beta_s * J * (J - 1.0) * FF_inv_trans;
        }
    }
    return;
}

void
MechanicsModel::compute_deviatoric_projection(TensorValue<double>& PP, const TensorValue<double>& FF)
{
    // compute deviatoric projection
    // PP_dev = PP - 1/3*(PP:FF)FF^-T
    PP = PP - (1.0/3.0)*PP.contract(FF)*FF.transpose().inverse();

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_neohookean(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& /*var_data*/,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const double J = FF.det();
    const double I1 = CC.tr();
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();

    const double mu_s = part_context->neohookean_parms.at(subdomain_ID)->mu_s;

    PP.zero();
    PP = mu_s * fast_pow_n23(J) * (FF - (1.0 / 3.0) * I1 * FF_inv_trans);

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_exp_neohookean(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& /*var_data*/,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const double J = FF.det();
    const double I1 = CC.tr();
    const double I1_bar = fast_pow_n23(J) * I1;
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int subdomain_ID = elem->subdomain_id();

    const double a = part_context->exp_neohookean_parms.at(subdomain_ID)->a;
    const double b = part_context->exp_neohookean_parms.at(subdomain_ID)->b;

    PP.zero();
    // the strain energy is W (I1) = (a/2b) * exp(b(I1 - 3))
    const double dW_dI1_bar = (0.5 * a) * exp(b * (I1_bar - 3.0));
    PP = dW_dI1_bar * dI1_bar_dFF(FF);

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_Guccione(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    TBOX_ASSERT(0 < var_data.size());
    TBOX_ASSERT(3 <= var_data[0]->size());
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();

    const double C = part_context->guccione_parms.at(subdomain_ID)->C;
    const double bf = part_context->guccione_parms.at(subdomain_ID)->bf;
    const double bt = part_context->guccione_parms.at(subdomain_ID)->bt;
    const double bfs = part_context->guccione_parms.at(subdomain_ID)->bfs;
    const unsigned int circ_fiber_index = part_context->guccione_parms.at(subdomain_ID)->circ_fiber_index;
    const bool DEV_projection = part_context->guccione_parms.at(subdomain_ID)->DEV_projection;

    VectorValue<double> f0;
    f0(0) = (*var_data[0])[0 + NDIM * circ_fiber_index];
    f0(1) = (*var_data[0])[1 + NDIM * circ_fiber_index];
    f0(2) = (*var_data[0])[2 + NDIM * circ_fiber_index];
    f0 = f0.unit();

    // stuff for logging min/max of the invariants
    const VectorValue<double> f = FF * f0;
    part_context->min_I4_one.at(subdomain_ID) = std::min(part_context->min_I4_one.at(subdomain_ID), f * f);
    part_context->max_I4_one.at(subdomain_ID) = std::max(part_context->max_I4_one.at(subdomain_ID), f * f);

    PP.zero();
    PP = Guccione_PK1_stress(bf, bt, bfs, C, f0, FF);

    if (DEV_projection)
    {
        compute_deviatoric_projection(PP, FF);
    }

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_HO(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double data_time,
    void* ctx)
{
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();

    TBOX_ASSERT(0 < var_data.size());
    if (part_context->ho_parms.at(subdomain_ID)->single_fiber_family)
    {
        TBOX_ASSERT(3 <= var_data[0]->size());
    }
    else
    {
        TBOX_ASSERT(6 <= var_data[0]->size());
    }
    const double mu_e = part_context->ho_parms.at(subdomain_ID)->mu_e;
    const double a = part_context->ho_parms.at(subdomain_ID)->a;
    const double b = part_context->ho_parms.at(subdomain_ID)->b;
    const double af = part_context->ho_parms.at(subdomain_ID)->af;
    const double bf = part_context->ho_parms.at(subdomain_ID)->bf;
    const double as = part_context->ho_parms.at(subdomain_ID)->as;
    const double bs = part_context->ho_parms.at(subdomain_ID)->bs;
    const double afs = part_context->ho_parms.at(subdomain_ID)->afs;
    const double bfs = part_context->ho_parms.at(subdomain_ID)->bfs;
    const double kf = part_context->ho_parms.at(subdomain_ID)->kf;
    const double ks = part_context->ho_parms.at(subdomain_ID)->ks;
    const bool enable_prestrain = part_context->ho_parms.at(subdomain_ID)->enable_prestrain;
    const double circ_prestrain_coeff = part_context->ho_parms.at(subdomain_ID)->circ_prestrain_coeff;
    const double radial_prestrain_coeff = part_context->ho_parms.at(subdomain_ID)->radial_prestrain_coeff;
    const unsigned int circ_fiber_index = part_context->ho_parms.at(subdomain_ID)->circ_fiber_index;
    const unsigned int radial_fiber_index = part_context->ho_parms.at(subdomain_ID)->radial_fiber_index;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->t_ramp;
    const bool use_scaled_invariants = part_context->ho_parms.at(subdomain_ID)->use_scaled_invariants;
    const bool turn_off_fibers_in_compression = part_context->ho_parms.at(subdomain_ID)->turn_off_fibers_in_compression;
    const bool DEV_projection = part_context->ho_parms.at(subdomain_ID)->DEV_projection;

    VectorValue<double> f0;
    f0(0) = (*var_data[0])[0 + NDIM * circ_fiber_index];
    f0(1) = (*var_data[0])[1 + NDIM * circ_fiber_index];
    f0(2) = (*var_data[0])[2 + NDIM * circ_fiber_index];
    f0 = f0.unit();

    // stuff for logging min/max of the invariants
    const VectorValue<double> f = FF * f0;
    part_context->min_I4_one.at(subdomain_ID) = std::min(part_context->min_I4_one.at(subdomain_ID), f * f);
    part_context->max_I4_one.at(subdomain_ID) = std::max(part_context->max_I4_one.at(subdomain_ID), f * f);

    VectorValue<double> s0(1.0, 0.0, 0.0);
    if (!part_context->ho_parms.at(subdomain_ID)->single_fiber_family)
    {
        s0(0) = (*var_data[0])[0 + NDIM * radial_fiber_index];
        s0(1) = (*var_data[0])[1 + NDIM * radial_fiber_index];
        s0(2) = (*var_data[0])[2 + NDIM * radial_fiber_index];
        s0 = s0.unit();

        // stuff for logging min/max of the invariants
        const VectorValue<double> s = FF * s0;
        part_context->min_I4_two.at(subdomain_ID) = std::min(part_context->min_I4_two.at(subdomain_ID), s * s);
        part_context->max_I4_two.at(subdomain_ID) = std::max(part_context->max_I4_two.at(subdomain_ID), s * s);
    }

    // define incompressible deformation in f0 direction
    // from reference configuration to initial configuration
    TensorValue<double> FFf0(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    TensorValue<double> FFs0(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    TensorValue<double> FF_total(FF);
    if (enable_prestrain)
    {
        TensorValue<double> eye(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        TensorValue<double> f0_prod_f0;
        outer_product(f0_prod_f0, f0, f0, NDIM);
        const double circ_prestrain =
        data_time < t_ramp ? (data_time * (circ_prestrain_coeff - 1.0) / t_ramp) : (circ_prestrain_coeff - 1.0);
        FFf0 = (1.0 + circ_prestrain) * f0_prod_f0 + (1 / sqrt(1.0 + circ_prestrain)) * (eye - f0_prod_f0);

        if (!part_context->ho_parms.at(subdomain_ID)->single_fiber_family)
        {
            // define incompressible deformation in s0 direction
            // from reference configuration to initial configuration
            TensorValue<double> s0_prod_s0;
            outer_product(s0_prod_s0, s0, s0, NDIM);
            const double radial_prestrain =
            data_time < t_ramp ? (data_time * (radial_prestrain_coeff - 1.0) / t_ramp) : (radial_prestrain_coeff - 1.0);
            FFs0 = (1.0 + radial_prestrain) * s0_prod_s0 + (1 / sqrt(1.0 + radial_prestrain)) * (eye - s0_prod_s0);
        }
        // incorporate prestrain into deformation gradient
        FF_total = FF * FFf0 * FFs0;
    }

    PP.zero();
    if (use_scaled_invariants)
    {
        PP = HO_PK1_stress_with_scaled_invariants(mu_e, a, b, af, bf, as, bs, afs, bfs, kf, ks, f0, s0, FF_total, turn_off_fibers_in_compression);
    }
    else
    {
        PP = HO_PK1_stress(mu_e, a, b, af, bf, as, bs, afs, bfs, kf, ks, f0, s0, FF_total, turn_off_fibers_in_compression);
    }

    // define the stress in the initial configuration (not the reference configuration)
    PP = PP * FFs0.transpose() * FFf0.transpose();

    if (DEV_projection)
    {
        compute_deviatoric_projection(PP, FF);
    }

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_biaxial_fung_type(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();

    TBOX_ASSERT(0 < var_data.size());
    TBOX_ASSERT(6 <= var_data[0]->size());

    const double C = part_context->biaxial_fung_type_parms.at(subdomain_ID)->C;
    const double A1 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A1;
    const double A2 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A2;
    const double A3 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A3;
    const double A4 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A4;
    const double A5 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A5;
    const double A6 = part_context->biaxial_fung_type_parms.at(subdomain_ID)->A6;
    const unsigned int circ_fiber_index = part_context->biaxial_fung_type_parms.at(subdomain_ID)->circ_fiber_index;
    const unsigned int radial_fiber_index = part_context->biaxial_fung_type_parms.at(subdomain_ID)->radial_fiber_index;
    const bool DEV_projection = part_context->biaxial_fung_type_parms.at(subdomain_ID)->DEV_projection;

    VectorValue<double> f0;
    f0(0) = (*var_data[0])[0 + NDIM * circ_fiber_index];
    f0(1) = (*var_data[0])[1 + NDIM * circ_fiber_index];
    f0(2) = (*var_data[0])[2 + NDIM * circ_fiber_index];
    f0 = f0.unit();

    VectorValue<double> s0;
    s0(0) = (*var_data[0])[0 + NDIM * radial_fiber_index];
    s0(1) = (*var_data[0])[1 + NDIM * radial_fiber_index];
    s0(2) = (*var_data[0])[2 + NDIM * radial_fiber_index];
    s0 = s0.unit();

    // stuff for logging min/max of the invariants
    const VectorValue<double> f = FF * f0;
    part_context->min_I4_one.at(subdomain_ID) = std::min(part_context->min_I4_one.at(subdomain_ID), f * f);
    part_context->max_I4_one.at(subdomain_ID) = std::max(part_context->max_I4_one.at(subdomain_ID), f * f);
    const VectorValue<double> s = FF * s0;
    part_context->min_I4_two.at(subdomain_ID) = std::min(part_context->min_I4_two.at(subdomain_ID), s * s);
    part_context->max_I4_two.at(subdomain_ID) = std::max(part_context->max_I4_two.at(subdomain_ID), s * s);

    PP.zero();
    PP = biaxial_fung_type_PK1_stress(C, A1, A2, A3, A4, A5, A6, f0, s0, FF);

    if (DEV_projection)
    {
        compute_deviatoric_projection(PP, FF);
    }

    return;
}

// https://doi.org/10.1007/s10237-019-01268-5
void
MechanicsModel::PK1_dev_stress_function_general_augustin(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();

    TBOX_ASSERT(0 < var_data.size());

    const double kappa = part_context->augustin_parms.at(subdomain_ID)->kappa;
    const double a = part_context->augustin_parms.at(subdomain_ID)->a;
    const double b = part_context->augustin_parms.at(subdomain_ID)->b;
    const double af = part_context->augustin_parms.at(subdomain_ID)->af;
    const double bf = part_context->augustin_parms.at(subdomain_ID)->bf;
    const unsigned int circ_fiber_index = part_context->augustin_parms.at(subdomain_ID)->circ_fiber_index;
    const bool DEV_projection = part_context->augustin_parms.at(subdomain_ID)->DEV_projection;

    TBOX_ASSERT(NDIM * (circ_fiber_index + 1) <= var_data[0]->size());

    VectorValue<double> f0;
    f0(0) = (*var_data[0])[0 + NDIM * circ_fiber_index];
    f0(1) = (*var_data[0])[1 + NDIM * circ_fiber_index];
    f0(2) = (*var_data[0])[2 + NDIM * circ_fiber_index];
    f0 = f0.unit();

    // stuff for logging min/max of the invariants
    const VectorValue<double> f = FF * f0;
    part_context->min_I4_one.at(subdomain_ID) = std::min(part_context->min_I4_one.at(subdomain_ID), f * f);
    part_context->max_I4_one.at(subdomain_ID) = std::max(part_context->max_I4_one.at(subdomain_ID), f * f);

    PP.zero();
    PP = augustin_PK1_stress(kappa, a, b, af, bf, f0, FF);

    if (DEV_projection)
    {
        compute_deviatoric_projection(PP, FF);
    }

    return;
}

void
MechanicsModel::PK1_dev_stress_function_general_nonlinear_spring(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*x*/, // current location
    const libMesh::Point& /*X*/, // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double data_time,
    void* ctx)
{
    PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const libMesh::subdomain_id_type subdomain_ID = elem->subdomain_id();
    const double af = part_context->nonlinear_spring_parms.at(subdomain_ID)->af;
    const double bf = part_context->nonlinear_spring_parms.at(subdomain_ID)->bf;
    const bool enable_prestrain = part_context->nonlinear_spring_parms.at(subdomain_ID)->enable_prestrain;
    const double prestrain_coeff = part_context->nonlinear_spring_parms.at(subdomain_ID)->prestrain_coeff;
    const unsigned int fiber_index = part_context->nonlinear_spring_parms.at(subdomain_ID)->fiber_index;
    const double t_ramp = part_context->general_parms.at(subdomain_ID)->t_ramp;
    const bool DEV_projection = part_context->nonlinear_spring_parms.at(subdomain_ID)->DEV_projection;

    TBOX_ASSERT(0 < var_data.size());
    TBOX_ASSERT(3 <= var_data[0]->size());

    VectorValue<double> f0;
    f0(0) = (*var_data[0])[0 + NDIM * fiber_index];
    f0(1) = (*var_data[0])[1 + NDIM * fiber_index];
    f0(2) = (*var_data[0])[2 + NDIM * fiber_index];
    f0 = f0.unit();

    // stuff for logging min/max of the invariants
    const VectorValue<double> f = FF * f0;
    part_context->min_I4_one.at(subdomain_ID) = std::min(part_context->min_I4_one.at(subdomain_ID), f * f);
    part_context->max_I4_one.at(subdomain_ID) = std::max(part_context->max_I4_one.at(subdomain_ID), f * f);

    // define incompressible deformation in f0 direction
    // from reference configuration to initial configuration
    TensorValue<double> FFf0(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    TensorValue<double> FF_total(FF);
    if (enable_prestrain)
    {
        TensorValue<double> eye(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        TensorValue<double> f0_prod_f0;
        outer_product(f0_prod_f0, f0, f0, NDIM);
        const double prestrain =
        data_time < t_ramp ? (data_time * (prestrain_coeff - 1.0) / t_ramp) : (prestrain_coeff - 1.0);
        FFf0 = (1.0 + prestrain) * f0_prod_f0 + (1 / sqrt(1.0 + prestrain)) * (eye - f0_prod_f0);

        // incorporate prestrain into deformation gradient
        FF_total = FF * FFf0;
    }

    PP.zero();
    PP = nonlinear_spring_PK1_stress(af, bf, f0, FF_total);

    // define the stress in the initial configuration (not the reference configuration)
    PP = PP * FFf0.transpose();

    if (DEV_projection)
    {
        compute_deviatoric_projection(PP, FF);
    }
    return;
}
