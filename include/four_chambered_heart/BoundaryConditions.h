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

#ifndef four_chambered_heart_boundary_conditions_h
#define four_chambered_heart_boundary_conditions_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// C++ includes
#include <vector>

// BoundaryConditions is a static class that provides data and functions
// required to implement loading conditions and boundary constraints for the LV
// model.
class BoundaryConditions
{
public:
    static double P_load, t_load;
    static double kappa_surface_tether;
    static double kappa_body_tether;
    static double kappa_surface_penalty;
    static double kappa_body_penalty;
    static double kappa_pericardium;
    static double eta_pericardium;

    static std::vector<BoundaryInfo*> boundary_info;

    static double loading_pressure(double time);

    static void
    loading_force_function(double& P,
                           const libMesh::VectorValue<double>& n,
                           const libMesh::VectorValue<double>& N,
                           const libMesh::TensorValue<double>& FF,
                           const libMesh::Point& x,
                           const libMesh::Point& X,
                           libMesh::Elem* elem,
                           unsigned short int side,
                           const std::vector<const std::vector<double>*>& system_var_data,
                           const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
                           double data_time,
                           void* ctx);

    static void tether_body_force_function(VectorValue<double>& F,
                                           const TensorValue<double>& FF,
                                           const libMesh::Point& x,
                                           const libMesh::Point& X,
                                           Elem* const elem,
                                           const vector<const vector<double>*>& var_data,
                                           const vector<const vector<VectorValue<double> >*>& grad_var_data,
                                           double time,
                                           void* ctx);

    static void tether_surface_force_function(
        VectorValue<double>& F,
        const libMesh::VectorValue<double>& n,
        const libMesh::VectorValue<double>& N,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        libMesh::Elem* elem,
        unsigned short int side,
        const std::vector<const std::vector<double>*>& system_var_data,
        const std::vector<const std::vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double data_time,
        void* ctx);

private:
    BoundaryConditions();
    BoundaryConditions(BoundaryConditions&);
    ~BoundaryConditions();
    BoundaryConditions& operator=(BoundaryConditions&);
};
#endif // four_chambered_heart_boundary_conditions_h
