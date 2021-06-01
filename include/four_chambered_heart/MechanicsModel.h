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

#ifndef four_chambered_heart_mechanics_model_h
#define four_chambered_heart_mechanics_model_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <four_chambered_heart/ModelParameters.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// helper function that is much faster than std::pow(a, -2.0/3.0). std::pow
// has trouble with some arguments which results in calling slowpow. Since we
// only call this function on values near 1 (Jacobians) we can be a bit more
// efficient and compute a*a first since we know that it won't overflow.
inline
double fast_pow_n23(const double a)
{
    return 1.0/std::cbrt(a*a);
}


// MechanicsModel is a static class that provides data and functions required to
// implement the heart mechanics model.
class MechanicsModel
{
public:
    static std::string HO_NAME;
    static std::string HGO_NAME;
    static std::string GUCCIONE_NAME;
    static std::string NEOHOOKEAN_NAME;
    static std::string EXP_NEOHOOKEAN_NAME;
    static std::string BIAXIAL_FUNG_TYPE_NAME;
    static std::string AUGUSTIN_NAME;
    static std::string NONLINEAR_SPRING_NAME;

    static int HO_ID;
    static int HGO_ID;
    static int GUCCIONE_ID;
    static int NEOHOOKEAN_ID;
    static int EXP_NEOHOOKEAN_ID;
    static int BIAXIAL_FUNG_TYPE_ID;
    static int AUGUSTIN_ID;
    static int NONLINEAR_SPRING_ID;

    static std::string ORIGINAL_ACTIVE_TENSION_FUNCTION_NAME;
    static std::string AUGUSTIN_ACTIVE_TENSION_FUNCTION_NAME;
    static std::string DATA_ACTIVE_TENSION_FUNCTION_NAME;
    static std::string PRESSURE_ACTIVE_TENSION_FUNCTION_NAME;

    static std::string LOG_VOLUMETRIC_ENERGY_NAME;
    static std::string QUADRATIC_VOLUMETRIC_ENERGY_NAME;

    static int LOG_VOLUMETRIC_ENERGY_ID;
    static int QUADRATIC_VOLUMETRIC_ENERGY_ID;

    static bool inverted_element;

    static void get_PK1_dev_stress_function_systems(vector<unsigned int>& systems);

    static double original_active_tension_function(const double time,
                                                   Elem* const elem,
                                                   void* ctx);

    static double augustin_active_tension_function(const double time,
                                                   Elem* const elem,
                                                   void* ctx);

    static double data_activation_function(const double time,
                                           Elem* const elem,
                                           void* ctx);

    static double pressure_activation_function(const double time,
                                               Elem* const elem,
                                               void* ctx);


    static void PK1_dev_stress_function_general_neohookean(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const std::vector<const std::vector<double>*>& var_data,
        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
        double data_time,
        void* ctx);

    static void PK1_dev_stress_function_general_exp_neohookean(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const std::vector<const std::vector<double>*>& var_data,
        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
        double data_time,
        void* ctx);

    static void
    PK1_dev_stress_function_general_HO(TensorValue<double>& PP,
                                       const TensorValue<double>& FF,
                                       const libMesh::Point& x,
                                       const libMesh::Point& X,
                                       Elem* const elem,
                                       const std::vector<const std::vector<double>*>& var_data,
                                       const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                       double data_time,
                                       void* ctx);

    static void
    PK1_dev_stress_function_general_Guccione(TensorValue<double>& PP,
                                             const TensorValue<double>& FF,
                                             const libMesh::Point& x,
                                             const libMesh::Point& X,
                                             Elem* const elem,
                                             const std::vector<const std::vector<double>*>& var_data,
                                             const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                             double data_time,
                                             void* ctx);

    static void PK1_dev_stress_function_general_biaxial_fung_type(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const std::vector<const std::vector<double>*>& var_data,
        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
        double data_time,
        void* ctx);

    static void
    PK1_dev_stress_function_general_augustin(TensorValue<double>& PP,
                                             const TensorValue<double>& FF,
                                             const libMesh::Point& x,
                                             const libMesh::Point& X,
                                             Elem* const elem,
                                             const std::vector<const std::vector<double>*>& var_data,
                                             const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                             double data_time,
                                             void* ctx);

    static void
    PK1_dev_stress_function_general_nonlinear_spring(TensorValue<double>& PP,
                                                     const TensorValue<double>& FF,
                                                     const libMesh::Point& x,
                                                     const libMesh::Point& X,
                                                     Elem* const elem,
                                                     const std::vector<const std::vector<double>*>& var_data,
                                                     const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                                     double data_time,
                                                     void* ctx);

    static void PK1_dev_stress_function(TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& x,
                                        const libMesh::Point& X,
                                        Elem* const elem,
                                        const std::vector<const std::vector<double>*>& var_data,
                                        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                        double data_time,
                                        void* ctx);

    static void PK1_dil_stress_function(TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& x,
                                        const libMesh::Point& X,
                                        Elem* const elem,
                                        const std::vector<const std::vector<double>*>& var_data,
                                        const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                        double data_time,
                                        void* ctx);

    static void PK1_combined_stress_function(TensorValue<double>& PP,
                                             const TensorValue<double>& FF,
                                             const libMesh::Point& x,
                                             const libMesh::Point& X,
                                             Elem* const elem,
                                             const std::vector<const std::vector<double>*>& var_data,
                                             const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
                                             double data_time,
                                             void* ctx);

    static void compute_deviatoric_projection(TensorValue<double>& PP, const TensorValue<double>& FF);

private:
    MechanicsModel();
    MechanicsModel(MechanicsModel&);
    ~MechanicsModel();
    MechanicsModel& operator=(MechanicsModel&);
};

// some helper functions for computing derivatives of invariants
// with respect to the deformation gradient

inline TensorValue<double>
dI1_dFF(const TensorValue<double>& FF)
{
    // I1 = I1(CC) = tr(FF^T FF)
    return 2.0 * FF;
}

inline TensorValue<double>
dI1_bar_dFF(const TensorValue<double>& FF)
{
    // I1_bar = I1(CC_bar) = tr(FF_bar^T FF_bar)
    // FF_bar = J^(-1/3) FF ===> det(FF_bar) = 1, I1_bar = J^(-2/3) I1
    const double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    return 2.0 * fast_pow_n23(J) * (FF - (1.0 / 3.0) * I1 * tensor_inverse_transpose(FF));
}

inline TensorValue<double>
dI4f_dFF(const TensorValue<double>& FF, const VectorValue<double>& f0)
{
    // I4f = f0 * CC * f0 = f0 * FF^T FF * f0 = (FF f0) * (FF f0)
    const VectorValue<double> f = FF * f0;
    return 2.0 * outer_product(f, f0);
}

inline TensorValue<double>
dI4f_bar_dFF(const TensorValue<double>& FF, const VectorValue<double>& f0)
{
    // I4f_bar = f0 * CC_bar * f0 = f0 * FF_bar^T FF_bar * f0 = (FF_bar f0) * (FF_bar f0)
    // FF_bar = J^(-1/3) FF ===> det(FF_bar) = 1, I4f_bar = J^(-2/3) I4f
    const double J = FF.det();
    const VectorValue<double> f = FF * f0;
    const double I4f = f * f;
    return 2.0 * fast_pow_n23(J) * (outer_product(f, f0) - (1.0 / 3.0) * I4f * tensor_inverse_transpose(FF));
}

inline TensorValue<double>
dI8fs_dFF(const TensorValue<double>& FF, const VectorValue<double>& f0, const VectorValue<double>& s0)
{
    // I8fs= f0 * CC * s0 = f0 * FF^T FF * s0 = (FF f0) * (FF s0)
    const VectorValue<double> f = FF * f0;
    const VectorValue<double> s = FF * s0;
    return outer_product(f, s0) + outer_product(s, f0);
}

inline TensorValue<double>
dI8fs_bar_dFF(const TensorValue<double>& FF, const VectorValue<double>& f0, const VectorValue<double>& s0)
{
    // I8fs_bar = f0 * CC_bar * s0 = f0 * FF_bar^T FF_bar * s0 = (FF_bar f0) * (FF_bar s0)
    // FF_bar = J^(-1/3) FF ===> det(FF_bar) = 1, I8fs_bar = J^(-2/3) I8fs
    const double J = FF.det();
    const VectorValue<double> f = FF * f0;
    const VectorValue<double> s = FF * s0;
    const double I8fs = f * s;
    return fast_pow_n23(J) *
           (outer_product(f, s0) + outer_product(s, f0) - (2.0 / 3.0) * I8fs * tensor_inverse_transpose(FF));
}

inline TensorValue<double>
dJ_dFF(const TensorValue<double>& FF)
{
    const double J = FF.det();
    return J * tensor_inverse_transpose(FF);
}

// helper functions for different material models
inline TensorValue<double>
HO_PK1_stress_with_scaled_invariants(const double mu_e,
                                     const double a,
                                     const double b,
                                     const double af,
                                     const double bf,
                                     const double as,
                                     const double bs,
                                     const double afs,
                                     const double bfs,
                                     const double kf,
                                     const double ks,
                                     const VectorValue<double> f0,
                                     const VectorValue<double> s0,
                                     const TensorValue<double> FF,
                                     const bool turn_off_fibers_in_compression)
{
    TensorValue<double> PP;
    PP.zero();

    const TensorValue<double> CC = FF.transpose() * FF;
    const double J = FF.det();
    const double J_n23 = fast_pow_n23(J);
    const double I1 = CC.tr();
    const double I1_bar = J_n23 * I1;
    const double I4f = (CC * f0) * f0;
    const double I4f_bar = J_n23 * I4f;
    const double I4s = (CC * s0) * s0;
    const double I4s_bar = J_n23 * I4s;
    const double I8fs = (CC * s0) * f0;
    const double I8fs_bar = J_n23 * I8fs;

    // compute I4_bar_star to account for fiber dispersion
    const double I4f_bar_temp = turn_off_fibers_in_compression ? std::max(I4f_bar, 1.0) : I4f_bar;
    const double I4s_bar_temp = turn_off_fibers_in_compression ? std::max(I4s_bar, 1.0) : I4s_bar;
    const double I4f_bar_star = kf * I1_bar + (1.0 - 3.0 * kf) * I4f_bar_temp;
    const double I4s_bar_star = ks * I1_bar + (1.0 - 3.0 * ks) * I4s_bar_temp;
    const TensorValue<double> dI4f_star_dFF = kf * dI1_bar_dFF(FF) + (1.0 - 3.0 * kf ) * dI4f_bar_dFF(FF, f0);
    const TensorValue<double> dI4s_star_dFF = ks * dI1_bar_dFF(FF) + (1.0 - 3.0 * ks ) * dI4f_bar_dFF(FF, s0);

    // the strain energy is
    // W (I1, I4i, I8fs) = (a/2b) * exp(b(I1 - 3))
    // + \sum_i (ai/2bi) ( exp( bi(I4i - 1)^2 ) - 1 )
    // + (afs/2bfs) ( exp( bfs I8fs^2 )  - 1 )

    const double dW_dI1_bar = (0.5 * a) * exp(b * (I1_bar - 3.0));
    const double dW_dI4f_bar_star = af * (I4f_bar_star - 1.0) * exp(bf * (I4f_bar_star - 1.0) * (I4f_bar_star - 1.0));
    const double dW_dI4s_bar_star = as * (I4s_bar_star - 1.0) * exp(bs * (I4s_bar_star - 1.0) * (I4s_bar_star - 1.0));
    const double dW_dI8fs_bar = afs * I8fs_bar * exp(bfs * I8fs_bar * I8fs_bar);

    PP += 2 * mu_e * (FF - tensor_inverse_transpose(FF));
    PP += dW_dI1_bar * dI1_bar_dFF(FF) + dW_dI4f_bar_star * dI4f_star_dFF + dW_dI4s_bar_star * dI4s_star_dFF +
          dW_dI8fs_bar * dI8fs_bar_dFF(FF, f0, s0);

    return PP;
}

inline TensorValue<double>
HO_PK1_stress(const double mu_e,
              const double a,
              const double b,
              const double af,
              const double bf,
              const double as,
              const double bs,
              const double afs,
              const double bfs,
              const double kf,
              const double ks,
              const VectorValue<double> f0,
              const VectorValue<double> s0,
              const TensorValue<double> FF,
              const bool turn_off_fibers_in_compression)
{
    TensorValue<double> PP;
    PP.zero();

    const TensorValue<double> CC = FF.transpose() * FF;
    const double I1 = CC.tr();
    const double I4f = (CC * f0) * f0;
    const double I4s = (CC * s0) * s0;
    const double I8fs = (CC * s0) * f0;

    // compute I4_bar_star to account for fiber dispersion
    const double I4f_temp = turn_off_fibers_in_compression ? std::max(I4f, 1.0) : I4f;
    const double I4s_temp = turn_off_fibers_in_compression ? std::max(I4s, 1.0) : I4s;
    const double I4f_star = kf * I1 + (1.0 - 3.0 * kf) * I4f_temp;
    const double I4s_star = ks * I1 + (1.0 - 3.0 * ks) * I4s_temp;
    const TensorValue<double> dI4f_star_dFF = kf * dI1_dFF(FF) + (1.0 - 3.0 * kf ) * dI4f_dFF(FF, f0);
    const TensorValue<double> dI4s_star_dFF = ks * dI1_dFF(FF) + (1.0 - 3.0 * ks ) * dI4f_dFF(FF, s0);

    // the strain energy is
    // W (I1, I4i, I8fs) = (a/2b) * exp(b(I1 - 3))
    // + \sum_i (ai/2bi) ( exp( bi(I4i - 1)^2 ) - 1 )
    // + (afs/2bfs) ( exp( bfs I8fs^2 )  - 1 )

    const double dW_dI1 = (0.5 * a) * exp(b * (I1- 3.0));
    const double dW_dI4f_star = af * (I4f_star - 1.0) * exp(bf * (I4f_star - 1.0) * (I4f_star - 1.0));
    const double dW_dI4s_star = as * (I4s_star - 1.0) * exp(bs * (I4s_star - 1.0) * (I4s_star - 1.0));
    const double dW_dI8fs = afs * I8fs* exp(bfs * I8fs * I8fs);

    PP += 2 * mu_e * (FF - tensor_inverse_transpose(FF));
    PP += dW_dI1 * dI1_dFF(FF) + dW_dI4f_star * dI4f_star_dFF + dW_dI4s_star * dI4s_star_dFF +
          dW_dI8fs * dI8fs_dFF(FF, f0, s0);

    return PP;
}


inline TensorValue<double>
Guccione_PK1_stress(const double bf,
                    const double bt,
                    const double bfs,
                    const double C,
                    const VectorValue<double> f0,
                    const TensorValue<double> FF0)
{
    TensorValue<double> PP;
    PP.zero();

    VectorValue<double> t00 = VectorValue<double>(-f0(2), 0, f0(0)).unit();
    VectorValue<double> t10 = f0.cross(t00);
    const TensorValue<double> RR(f0(0), f0(1), f0(2), t00(0), t00(1), t00(2), t10(0), t10(1), t10(2));
    const TensorValue<double> RR_t = RR.transpose();

    const TensorValue<double> FF = RR * FF0 * RR_t;
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const TensorValue<double> EE = 0.5 * (CC - II);

    // Passive contribution
    double Q = bf * EE(0, 0) * EE(0, 0) +
               bt * (EE(1, 1) * EE(1, 1) + EE(2, 2) * EE(2, 2) + EE(1, 2) * EE(1, 2) + EE(2, 1) * EE(2, 1)) +
               bfs * (EE(0, 1) * EE(0, 1) + EE(1, 0) * EE(1, 0) + EE(0, 2) * EE(0, 2) + EE(2, 0) * EE(2, 0));
    double CeQ = C * exp(Q);

    TensorValue<double> SS;
    SS(0, 0) = CeQ * bf * EE(0, 0);
    SS(0, 1) = CeQ * bfs * 0.5 * (EE(0, 1) + EE(1, 0));
    SS(0, 2) = CeQ * bfs * 0.5 * (EE(0, 2) + EE(2, 0));
    SS(1, 0) = CeQ * bfs * 0.5 * (EE(0, 1) + EE(1, 0));
    SS(1, 1) = CeQ * bt * EE(1, 1);
    SS(1, 2) = CeQ * bt * 0.5 * (EE(1, 2) + EE(2, 1));
    SS(2, 0) = CeQ * bfs * 0.5 * (EE(0, 2) + EE(2, 0));
    SS(2, 1) = CeQ * bt * 0.5 * (EE(1, 2) + EE(2, 1));
    SS(2, 2) = CeQ * bt * EE(2, 2);

    // PK-1 stress is computed in the physical coordinate system.
    const TensorValue<double> SS0 = RR_t * SS * RR;
    PP = FF0 * SS0;
    return PP;
}

inline TensorValue<double>
biaxial_fung_type_PK1_stress(const double C,
                             const double A1,
                             const double A2,
                             const double A3,
                             const double A4,
                             const double A5,
                             const double A6,
                             const VectorValue<double> f0,
                             const VectorValue<double> s0,
                             const TensorValue<double> FF0)
{
    TensorValue<double> PP;
    PP.zero();

    // generate the third orthonomal basis vector
    const VectorValue<double> fxs = f0.cross(s0);
    const TensorValue<double> RR(f0(0), f0(1), f0(2), s0(0), s0(1), s0(2), fxs(0), fxs(1), fxs(2));
    const TensorValue<double> RR_t = RR.transpose();

    // transform FF0
    const TensorValue<double> FF = RR * FF0 * RR_t;
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.1);
    const TensorValue<double> EE = 0.5 * (CC - II);

    const double Q = A1 * EE(0, 0) * EE(0, 0) + A2 * EE(1, 1) * EE(1, 1) + 2 * A3 * EE(0, 0) * EE(1, 1) +
                     A4 * EE(0, 1) * EE(0, 1) + 2 * A5 * EE(0, 0) * EE(0, 1) + 2 * A6 * EE(1, 1) * EE(0, 1);

    // generate stress tensor
    TensorValue<double> SS;
    SS(0, 0) = C * exp(Q) * (A1 * EE(0, 0) + A3 * EE(1, 1) + A5 * EE(0, 1));
    SS(1, 0) = C * exp(Q) * (A4 * EE(0, 1) + A5 * EE(0, 0) + A6 * EE(1, 1));
    SS(0, 1) = SS(1, 0);
    SS(1, 1) = C * exp(Q) * (A2 * EE(1, 1) + A3 * EE(0, 0) + A6 * EE(0, 1));

    // change back to original orientation
    const TensorValue<double> SS0 = RR_t * SS * RR;
    PP = FF0 * SS0;

    return PP;
}

inline TensorValue<double>
augustin_PK1_stress(const double kappa,
                    const double a,
                    const double b,
                    const double af,
                    const double bf,
                    const VectorValue<double> f0,
                    const TensorValue<double> FF)
{
    TensorValue<double> PP;
    PP.zero();

    const TensorValue<double> CC = FF.transpose() * FF;
    //const double J = FF.tr();
    const double I1 = CC.tr();
    const double I4 = (CC * f0) * f0;

    // W(I1, I4) = a/(2*b)*(exp(b(I1-3))-1) +
    // af/(2*bf)*(exp(bf(kappa*I1+(1-3*kappa)I4-1)^2)-1)

    PP += a * 0.5 * exp(b * (I1 - 3)) * dI1_dFF(FF);
    PP += af * exp(bf * pow(kappa * I1 + (1 - 3 * kappa) * I4 - 1, 2.0)) * (kappa * I1 + (1 - 3 * kappa) * I4 - 1) *
          (kappa * dI1_dFF(FF) + (1 - 3 * kappa) * dI4f_dFF(FF, f0));

    return PP;
}

inline TensorValue<double>
nonlinear_spring_PK1_stress(const double af,
                            const double bf,
                            const VectorValue<double> f0,
                            const TensorValue<double> FF)
{
    TensorValue<double> PP;
    PP.zero();

    const TensorValue<double> CC = FF.transpose() * FF;
    const double I4fm1 = (CC * f0) * f0 - 1.0;
    //const double lambda = std::sqrt(I4f);

    // the strain energy is
    // W (I4f) = 2 * af / ( bf+1 ) * ( sqrt(I4f) - 1 )^(bf+1)
    // note that this model assumes no compressive forces, i.e. PP = zero tensor when lambda <= 1.0
    //const double dW_dI4f = (lambda > 1.0) ? af * std::pow(lambda - 1.0, bf) / lambda : 0.0;
    //PP += dW_dI4f* dI4f_dFF(FF, f0);

    // W(I1, I4f) = af / 2 * (I1 - 3) + bf / 3 * (I4f - 1)^3 * (I4f>1)
    const double dW_dI1 = af / 2.0;
    const double dW_dI4f = (I4fm1 > 0) ? bf * I4fm1 * I4fm1 : 0.0;
    PP += dW_dI1 * dI1_dFF(FF) + dW_dI4f * dI4f_dFF(FF,f0);
    return PP;
}

#endif // four_chambered_heart_mechanics_model_h
