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

// LIBMESH INCLUDES
#include <libmesh/boundary_info.h>

// APPLICATION INCLUDES
#include <four_chambered_heart/BoundaryConditions.h>
#include <four_chambered_heart/CirculationModel.h>
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/ModelInitialization.h>
#include <four_chambered_heart/PartContext.h>
#include <four_chambered_heart/ModelParameters.h>

// STATIC VARIABLES
double BoundaryConditions::P_load, BoundaryConditions::t_load;
double BoundaryConditions::kappa_surface_tether;
double BoundaryConditions::kappa_body_tether;
double BoundaryConditions::kappa_surface_penalty;
double BoundaryConditions::kappa_body_penalty;
double BoundaryConditions::kappa_pericardium;
double BoundaryConditions::eta_pericardium;

std::vector<BoundaryInfo*> BoundaryConditions::boundary_info;

// CLASS IMPLEMENTATION

double
BoundaryConditions::loading_pressure(double time)
{
    return P_load * min(time, t_load) / t_load;
}

void
BoundaryConditions::loading_force_function(
    double& P,
    const libMesh::VectorValue<double>& /*n*/,
    const libMesh::VectorValue<double>& /*N*/,
    const libMesh::TensorValue<double>& /*FF*/,
    const libMesh::Point& /*x*/,
    const libMesh::Point& /*X*/,
    libMesh::Elem* elem,
    unsigned short int side,
    const std::vector<const std::vector<double>*>& /*system_var_data*/,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double data_time,
    void* ctx)
{
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int part = part_context->part_id;

    std::vector<boundary_id_type> bdry_ids;
    boundary_info[part]->boundary_ids(elem, side, bdry_ids);
    if (find_first_of(bdry_ids.begin(),
                      bdry_ids.end(),
                      MeshInfo::surface_pressure_IDs[part].begin(),
                      MeshInfo::surface_pressure_IDs[part].end()) != bdry_ids.end())
    {
        P = loading_pressure(data_time);
    }
    else
    {
        P = 0.0;
    }
    return;
}

void
BoundaryConditions::tether_body_force_function(VectorValue<double>& F,
                                               const TensorValue<double>& /*FF*/,
                                               const libMesh::Point& x,
                                               const libMesh::Point& X,
                                               Elem* const elem /*elem*/,
                                               const vector<const vector<double>*>& /*var_data*/,
                                               const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                               double /*time*/,
                                               void* ctx)
{
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int part = part_context->part_id;

    F.zero();
    const short int block_id = elem->subdomain_id();
    
    if (MeshInfo::circ_model_body_tethering)
    {
        if (find(MeshInfo::circ_model_body_force_IDs[part].begin(), MeshInfo::circ_model_body_force_IDs[part].end(), block_id) !=
                MeshInfo::circ_model_body_force_IDs[part].end())
        {
            const Afterload_Parms afterload_parms = MeshInfo::subdomain_to_afterload_parms[block_id];
            const int axis = afterload_parms.axis;
            const int side = afterload_parms.side;
            // the code below assumes the model is inside a square domain [-L, L]^3
            const bool apply_force = (side == 0) ? 
                (fabs(X(axis) + MeshInfo::Cartesian_L[side]) < 0.1 * MeshInfo::Cartesian_L[side]) : 
                 (fabs(X(axis) - MeshInfo::Cartesian_L[side]) < 0.1 * MeshInfo::Cartesian_L[side]);
            if (apply_force)
            {
                F = kappa_body_penalty * (X - x);
            }
        }
    }
    
    if (find(MeshInfo::tether_body_force_IDs[part].begin(), MeshInfo::tether_body_force_IDs[part].end(), block_id) !=
        MeshInfo::tether_body_force_IDs[part].end())
    {
        F = kappa_body_tether * (X - x);
    }

    return;
}

void
BoundaryConditions::tether_surface_force_function(
    VectorValue<double>& F,
    const libMesh::VectorValue<double>& n,
    const libMesh::VectorValue<double>& /*N*/,
    const libMesh::TensorValue<double>& /*FF*/,
    const libMesh::Point& x,
    const libMesh::Point& X,
    libMesh::Elem* elem,
    unsigned short int side,
    const std::vector<const std::vector<double>*>& system_var_data,
    const std::vector<const std::vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double /*data_time*/,
    void* ctx)
{
    const PartContext* part_context = reinterpret_cast<PartContext*>(ctx);
    const int part = part_context->part_id;

    F.zero();
    std::vector<boundary_id_type> bdry_ids;
    boundary_info[part]->boundary_ids(elem, side, bdry_ids);

    // general tethering force
    if(MeshInfo::circ_model_surface_tethering)
    {
        if (find_first_of(bdry_ids.begin(),
                bdry_ids.end(),
                MeshInfo::tether_surface_force_IDs[part].begin(),
                MeshInfo::tether_surface_force_IDs[part].end()) != bdry_ids.end())
        {
            F = kappa_surface_penalty * (X - x);
        }
    }

    // pericardial tethering force
    if (find_first_of(bdry_ids.begin(),
                      bdry_ids.end(),
                      MeshInfo::pericardial_tethering_IDs.begin(),
                      MeshInfo::pericardial_tethering_IDs.end()) != bdry_ids.end())
    {
        VectorValue<double> velocity((*system_var_data[0])[0], (*system_var_data[0])[1], (*system_var_data[0])[2]);
        F = n * (kappa_pericardium * (X - x) * n - eta_pericardium * velocity * n);
    }

    return;
}
