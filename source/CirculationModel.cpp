// Filename: CirculationModel.cpp
// Created on 20 Aug 2007 by Boyce Griffith

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>
#include <libmesh/point.h>

// C++ STDLIB INCLUDES
#include <cassert>

#include <Eigen/Dense>

#include <four_chambered_heart/CirculationModel.h>
#include <four_chambered_heart/MeshInfo.h>
#include <four_chambered_heart/ModelParameters.h>

using namespace Eigen;

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
    
// Conversion factors.
static const double flconv = 1 / 0.06;  // l/min ===> ml/sec
static const double prconv = 1333.2239; // mmHg  ===> dyne/cm^2

vector<double> data_time, pressure;
inline double
interpolate_data(const double& t, const vector<double>& data, const vector<double>& time)
{
    const double t_end = time.back();
    double t_mod = fmod(t, t_end);
    vector<double>::const_iterator posn = std::lower_bound(time.begin(), time.end(), t_mod);
    size_t k = posn - time.begin();
    if (t_mod < time[k]) k -= 1;
    const double t0 = time[k];
    const double t1 = time[k + 1];
    if (k >= time.size()) TBOX_ERROR("indexing error!\n");
    if (t_mod < t0) TBOX_ERROR("t0 error!\n");
    if (t_mod > t1) TBOX_ERROR("t1 error!\n");
    const double dt = t1 - t0;
    const double del = t_mod - t0;
    const double value0 = data[k];
    const double value1 = data[k + 1];
    return value0 + (value1 - value0) * del / dt;
} // interpolate_data
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel::CirculationModel(const string& object_name,
                                   Pointer<Database> input_db,
                                   const std::vector<Afterload_Parms> afterload_parms,
                                   bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_time(0.0),
      d_bdry_interface_level_number(numeric_limits<int>::max())
{
#if !defined(NDEBUG)
    assert(!object_name.empty());
#endif
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }
    if (input_db)
    {
        d_bdry_interface_level_number =
            input_db->getIntegerWithDefault("bdry_interface_level_number", d_bdry_interface_level_number);
        d_data_time_filename = input_db->getString("data_time_filename");
        
        // make directory for data files
        d_data_directory_name = input_db->getString("data_directory_name");
        Utilities::recursiveMkdir(d_data_directory_name);
        
        d_nsrc = afterload_parms.size();
        // resize vectors
        d_R_P.resize(d_nsrc);
        d_R_C.resize(d_nsrc);
        d_C.resize(d_nsrc);
        d_Q_sink.resize(d_nsrc);
        d_P_ramp.resize(d_nsrc);
        d_t_ramp.resize(d_nsrc);
        d_rsrc.resize(d_nsrc);
        d_q_prescribed_size.resize(d_nsrc);
        d_q_prescribed_times.resize(d_nsrc);
        d_q_prescribed_values.resize(d_nsrc);        
        d_use_velocity_bcs.resize(d_nsrc);
        d_use_Q_sink.resize(d_nsrc);
        d_src_axis.resize(d_nsrc);
        d_src_side.resize(d_nsrc);
        d_posn.resize(NDIM * d_nsrc);
        d_srcname.resize(d_nsrc);
        
        // bring in boolean determining if we want to test the circ models
        d_test_circ_models = input_db->getBoolWithDefault("test_circ_models", false);
        // bring in flow waveform if we want to test circ models
        if (d_test_circ_models)
        {
            d_input_flow_data_size = input_db->getInteger("input_flow_data_size");
            d_input_flow_times.resize(d_input_flow_data_size);
            d_input_flow_values.resize(d_input_flow_data_size);
            input_db->getDoubleArray("input_flow_times", &d_input_flow_times[0], d_input_flow_data_size);
            TBOX_ASSERT(!d_input_flow_times.empty());
            input_db->getDoubleArray("input_flow_values", &d_input_flow_values[0], d_input_flow_data_size);
            TBOX_ASSERT(!d_input_flow_values.empty());
        }
        
        // loop over afterload sources to bring in parameters and filenames
        for (int src_id = 0; src_id < d_nsrc; ++src_id)
        {
            d_srcname[src_id] = afterload_parms[src_id].nodeset_name;
            // we do not want to use velocity bcs if we are just testing the circ models
            if (d_test_circ_models && d_use_velocity_bcs[src_id])
            {
                TBOX_ERROR("CirculationModel::CirculationModel: velocity boundary conditions should"
                        " not be used for source " << d_srcname[src_id] << " since we are just testing"
                        " the circ models");
            }
 
            Pointer<Database> src_db = input_db->getDatabase(d_srcname[src_id]);
            d_rsrc[src_id] = afterload_parms[src_id].radius;
            d_use_velocity_bcs[src_id] = src_db->getBool("use_velocity_bcs");
            d_use_Q_sink[src_id] = src_db->getBool("use_Q_sink");
            if (d_use_velocity_bcs[src_id] && d_use_Q_sink[src_id])
            {
                TBOX_ERROR("CirculationModel::CirculationModel: cannot use both velocity"
                        " boundary conditions and Q_sink for afterload model with name: "
                        << d_srcname[src_id] << ".\n");
            }
            if (d_use_velocity_bcs[src_id])
            {
                d_q_prescribed_size[src_id] = src_db->getInteger("q_prescribed_size");
                d_q_prescribed_times[src_id].resize(d_q_prescribed_size[src_id]);
                d_q_prescribed_values[src_id].resize(d_q_prescribed_size[src_id]);
                src_db->getDoubleArray("q_prescribed_times", &d_q_prescribed_times[src_id][0], d_q_prescribed_size[src_id]);
                TBOX_ASSERT(!d_q_prescribed_times[src_id].empty());
                src_db->getDoubleArray("q_prescribed_values", &d_q_prescribed_values[src_id][0], d_q_prescribed_size[src_id]);
                TBOX_ASSERT(!d_q_prescribed_values[src_id].empty());
            }
            if (d_use_Q_sink[src_id])
            {
                d_Q_sink[src_id] = src_db->getDouble("Q_sink");
            }
            d_src_axis[src_id] = afterload_parms[src_id].axis;
            d_src_side[src_id] = afterload_parms[src_id].side;
            d_posn[src_id * NDIM] = afterload_parms[src_id].centroid(0);
            d_posn[src_id * NDIM + 1] = afterload_parms[src_id].centroid(1);
            d_posn[src_id * NDIM + 2] = afterload_parms[src_id].centroid(2);
            d_R_P[src_id] = src_db->getDouble("R_P");
            d_R_C[src_id] = src_db->getDouble("R_C");
            d_C[src_id] = src_db->getDouble("C");
            d_P_ramp[src_id] = src_db->getDouble("P_ramp");
            d_t_ramp[src_id] = src_db->getDouble("t_ramp");
        }

        // setup vector for location index to a vector of source ids
        d_location_index_to_src_id.resize(2 * NDIM);
        for (int src_id = 0; src_id < d_nsrc; ++src_id)
        {
            const int location_index = d_src_axis[src_id] * 2 + d_src_side[src_id];
            d_location_index_to_src_id[location_index].push_back(src_id);
        }
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    else
    {
        // resize vectors and initialize
        d_qsrc.resize(d_nsrc);
        d_psrc.resize(d_nsrc);
        d_usrc_prescribed.resize(d_nsrc);
        d_P_Wk.resize(d_nsrc);
        for (int src_id = 0; src_id < d_nsrc; ++src_id)
        {
            d_qsrc[src_id] = 0.0;
            d_psrc[src_id] = 0.0;
            d_usrc_prescribed[src_id] = 0.0;
            d_P_Wk[src_id] = 0.0;
        }
    }
    return;
} // CirculationModel

// Backward Euler update for windkessel model.
void
CirculationModel::windkessel_be_update(double& P, const int& src_id, const double& Q, const double& t, const double& dt)
{
    if (t < d_t_ramp[src_id])
    {
        if (d_use_Q_sink[src_id])
        {
            // the model that uses a Q_sink is the same as the one with Q_sink set to zero, 
            // except in this case R_p is set to infinity.
            d_P_Wk[src_id] = d_P_Wk[src_id] + (dt / d_C[src_id]) * (Q - (t / d_t_ramp[src_id]) * d_Q_sink[src_id]);
        }
        else
        {
            d_P_Wk[src_id] = (t + dt) * d_P_ramp[src_id] / d_t_ramp[src_id];
        }
    }
    else
    {
        if (d_use_Q_sink[src_id])
        {
            d_P_Wk[src_id] = d_P_Wk[src_id] + (dt / d_C[src_id]) * (Q - d_Q_sink[src_id]);
        }
        else
        {
            d_P_Wk[src_id] = (d_C[src_id] * d_P_Wk[src_id] + Q * dt) * d_R_P[src_id] / (d_C[src_id] * d_R_P[src_id] + dt);
        }
    }
    P = d_P_Wk[src_id] + d_R_C[src_id] * Q;
    return;
}

CirculationModel::~CirculationModel()
{
    return;
} // ~CirculationModel

void
CirculationModel::advanceTimeDependentData(const double dt,
                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int U_idx,
                                           const int /*P_idx*/,
                                           const int /*wgt_cc_idx*/,
                                           const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    std::fill(d_qsrc.begin(), d_qsrc.end(), 0.0);
    
    // actually sample the flow if we are NOT testing the circ models
    if(!d_test_circ_models) 
    {
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                if (pgeom->getTouchesRegularBoundary())
                {
                    Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                    Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                    const Box<NDIM>& patch_box = patch->getBox();
                    const double* const x_lower = pgeom->getXLower();
                    const double* const dx = pgeom->getDx();
                    double dV = 1.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        dV *= dx[d];
                    }
                    double X[NDIM];
                    
                    for (int ss = 0; ss < d_nsrc; ++ss)
                    {
                        const int axis = d_src_axis[ss];
                        const int side = d_src_side[ss];
                        if (pgeom->getTouchesRegularBoundary(axis, side))
                        {
                            const double rsrc = d_rsrc[ss];
                            IBTK::Point posn;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                posn(d) = d_posn[ss * NDIM + d];
                            }
                            
                            // normal vector
                            Vector n;
                            for (int nn = 0; nn < NDIM; ++nn)
                            {
                                n[nn] = (nn == axis ? pow(-1.0, (double)side + 1.0) : 0.0);
                            }
                            
                            Box<NDIM> side_box = patch_box;
                            if (side == 0)
                            {
                                side_box.lower(axis) = patch_box.lower(axis) + 1;
                                side_box.upper(axis) = patch_box.lower(axis) + 1;
                            }
                            else
                            {
                                side_box.lower(axis) = patch_box.upper(axis) + 1;
                                side_box.upper(axis) = patch_box.upper(axis) + 1;
                            }
                            
                            for (Box<NDIM>::Iterator b(side_box); b; b++)
                            {
                                const hier::Index<NDIM>& i = b();
                                double r_sq = 0.0;
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    X[d] =
                                            x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                    if (d != axis) r_sq += pow(X[d] - posn[d], 2.0);
                                }
                                const double r = sqrt(r_sq);
                                if (r <= rsrc)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = n[axis] * dV / dx[axis];
                                        d_qsrc[ss] += (*U_data)(i_s)*dA;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(&d_qsrc[0], d_nsrc);
    }

    const double t = d_time;
    // update pressures in all of the sources.
    for (int src_id = 0; src_id < d_nsrc; ++src_id)
    {
        if (d_test_circ_models) d_qsrc[src_id] = flow_function(t);
        double P;
        const double Q = d_qsrc[src_id]; // flow passed into function should be in mL/s
        // setting pressure boundary condition from windkessel model
        if (!d_use_velocity_bcs[src_id]) 
        {
            windkessel_be_update(P, src_id, Q, t, dt);
            d_psrc[src_id] = P * prconv; // pressure boundary condition should be in dynes/cm^2
        }
        else // setting velocity boundary condition from a specified flow
        {
            const double flow_value = interpolate_data(t, d_q_prescribed_values[src_id], d_q_prescribed_times[src_id]);
            d_usrc_prescribed[src_id] = flow_value / (M_PI * d_rsrc[src_id] * d_rsrc[src_id]);
        }
    }

    // Update the current time.
    d_time += dt;

    // Output the updated values.
    const long precision = plog.precision();
    plog.unsetf(ios_base::showpos);
    plog.unsetf(ios_base::scientific);

    plog << "============================================================================\n"
         << "Circulation model variables at time " << d_time << ":\n";

    for (int src_id = 0; src_id < d_nsrc; ++src_id)
    {
        plog << "afterload model name: " << d_srcname[src_id] << ": ";
        plog << "\n";
        plog << "Q = ";
        plog.setf(ios_base::showpos);
        plog.setf(ios_base::scientific);
        plog.precision(12);
        plog << d_qsrc[src_id] / flconv << "    ";

        if (!d_use_velocity_bcs[src_id])
        {
            plog << "P = ";
            plog.setf(ios_base::showpos);
            plog.setf(ios_base::scientific);
            plog.precision(12);
            plog << d_psrc[src_id] / prconv << "    ";
            
            plog << "P_Wk = ";
            plog.setf(ios_base::showpos);
            plog.setf(ios_base::scientific);
            plog.precision(12);
            plog << d_P_Wk[src_id] << "    ";
        }
        plog << "\n";
    }

    plog << "flow units: liter/min    ";
    plog << "pressure units: mmHg\n";

    plog << "============================================================================\n";

    plog.unsetf(ios_base::showpos);
    plog.unsetf(ios_base::scientific);
    plog.precision(precision);

    // Write the current state to disk.
    writeDataFile();
    return;
} // advanceTimeDependentData

void
CirculationModel::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_time", d_time);
    db->putInteger("d_nsrc", d_nsrc);
    db->putDoubleArray("d_qsrc", &d_qsrc[0], d_nsrc);
    db->putDoubleArray("d_psrc", &d_psrc[0], d_nsrc);
    db->putDoubleArray("d_usrc_prescribed", &d_usrc_prescribed[0], d_nsrc);
    db->putDoubleArray("d_P_Wk", &d_P_Wk[0], d_nsrc);
    db->putInteger("d_bdry_interface_level_number", d_bdry_interface_level_number);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CirculationModel::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool files_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        const unsigned int io_precision = 8;
        if (!files_initialized && !from_restart)
        {
            ofstream time_stream(d_data_time_filename.c_str(), ios::out);
            time_stream.precision(io_precision);
            time_stream << d_time << "\n";
            for (int src_id = 0; src_id < d_nsrc; ++src_id)
            {
                ofstream P_stream((d_data_directory_name + "/P_" + d_srcname[src_id] + ".dat").c_str(), ios::out);
                P_stream.precision(io_precision);
                ofstream PWk_stream((d_data_directory_name + "/PWk_" + d_srcname[src_id] + ".dat").c_str(), ios::out);
                PWk_stream.precision(io_precision);
                ofstream Q_stream((d_data_directory_name + "/Q_" + d_srcname[src_id] + ".dat").c_str(), ios::out);
                Q_stream.precision(io_precision);
                P_stream << d_psrc[src_id] / prconv << "\n";
                PWk_stream << d_P_Wk[src_id] << "\n";
                Q_stream << d_qsrc[src_id] << "\n";
            }
            files_initialized = true;
        }
        ofstream time_stream(d_data_time_filename.c_str(), ios::app);
        time_stream.precision(io_precision);
        time_stream << d_time << "\n";
        for (int src_id = 0; src_id < d_nsrc; ++src_id)
        {
            ofstream P_stream((d_data_directory_name + "/P_" + d_srcname[src_id] + ".dat").c_str(), ios::app);
            P_stream.precision(io_precision);
            ofstream PWk_stream((d_data_directory_name + "/PWk_" + d_srcname[src_id] + ".dat").c_str(), ios::app);
            PWk_stream.precision(io_precision);
            ofstream Q_stream((d_data_directory_name + "/Q_" + d_srcname[src_id] + ".dat").c_str(), ios::app);
            Q_stream.precision(io_precision);
            P_stream << d_psrc[src_id] / prconv << "\n";
            PWk_stream << d_P_Wk[src_id] << "\n";
            Q_stream << d_qsrc[src_id] << "\n";
        }
    }
    return;
} // writeDataFile

void
CirculationModel::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_time = db->getDouble("d_time");
    d_nsrc = db->getInteger("d_nsrc");
    d_qsrc.resize(d_nsrc);
    d_psrc.resize(d_nsrc);
    d_usrc_prescribed.resize(d_nsrc);
    d_P_Wk.resize(d_nsrc);
    if (db->keyExists("d_qsrc"))
    {
        db->getDoubleArray("d_qsrc", &d_qsrc[0], d_nsrc);
        TBOX_ASSERT(!d_qsrc.empty());
    }
    if(db->keyExists("d_psrc"))
    {
        db->getDoubleArray("d_psrc", &d_psrc[0], d_nsrc);
        TBOX_ASSERT(!d_psrc.empty());
    }
    if (db->keyExists("d_usrc_prescribed"))
    {
        db->getDoubleArray("d_usrc_prescribed", &d_usrc_prescribed[0], d_nsrc);
        TBOX_ASSERT(!d_usrc_prescribed.empty());
    }
    if(db->keyExists("d_P_Wk"))
    {
        db->getDoubleArray("d_P_Wk", &d_P_Wk[0], d_nsrc);
        TBOX_ASSERT(!d_P_Wk.empty());
    }
    d_bdry_interface_level_number = db->getInteger("d_bdry_interface_level_number");
    return;
} // getFromRestart

double 
CirculationModel::flow_function(double time)
{
    const double final_time = d_input_flow_times[d_input_flow_data_size-1];
    const double time_mod = fmod(time, final_time);
    // get iterator to the first data point larger than the current scaled time
    auto it_t1=std::lower_bound(d_input_flow_times.begin(), d_input_flow_times.end(), time_mod); //
    // get iterator to the previous data point
    auto it_t0=it_t1-1;

    // find the index in the array corresponding to the iterators
    int i1 = std::distance(d_input_flow_times.begin(), it_t1);
    int i0 = i1-1;

    // Get the time interval in which the scaled_t falls
    double t0 = *it_t0;
    double t1 = *it_t1;

    // get the (normalized) force interval
    double f0 = d_input_flow_values[i0];
    double f1 = d_input_flow_values[i1];

    // interpolate linearly
    double f =  (time_mod - t0) * (f1 - f0) / (t1 - t0) + f0;

    return f;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
