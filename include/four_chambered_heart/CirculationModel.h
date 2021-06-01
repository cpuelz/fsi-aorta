// Filename: CirculationModel.h
// Created on 04 May 2007 by Boyce Griffith

#ifndef four_chambered_heart_circulation_model_h
#define four_chambered_heart_circulation_model_h

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// IBTK INCLUDES
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <tbox/Database.h>
#include <tbox/Serializable.h>

// libMesh includes
#include <libmesh/point.h>

// C++ STDLIB INCLUDES
#include <vector>

// NAMESPACE
#include <ibamr/app_namespaces.h>

#include <four_chambered_heart/MeshInfo.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel
 */
class CirculationModel : public Serializable
{
public:
    /*!
     * \brief The object name.
     */
    string d_object_name;
    
    /*!
     * \brief Whether the object is registered with the restart manager.
     */
    bool d_registered_for_restart;
    
    /*!
     * \brief true if we want to test the circ models
     */
    bool d_test_circ_models;
    
    /*!
     * \brief number of data points we have for the input flow waveform
     */
    unsigned int d_input_flow_data_size; 
    
    /*!
     * \brief vectors for storing input flow waveform info if we are testing the 
     * circ models
     */
    vector<double> d_input_flow_times;
    vector<double> d_input_flow_values;

    /*!
     * \brief Windkessel model data.
     */
    double d_time;
    int d_nsrc;
    vector<double> d_qsrc, d_psrc, d_usrc_prescribed, d_rsrc, d_P_Wk;
    
    // vectors for bringing in flow waveforms
    vector<unsigned int> d_q_prescribed_size;
    vector<vector<double> > d_q_prescribed_times, d_q_prescribed_values;
    
    vector<bool> d_use_velocity_bcs;
    vector<bool> d_use_Q_sink;
    vector<int> d_src_axis;
    vector<int> d_src_side;
    vector<vector<int> > d_location_index_to_src_id;
    vector<double> d_posn;
    vector<string> d_srcname;

    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief filename containing the time data
     */
    std::string d_data_time_filename;
    
    /*!
     * \brief string containing the directory for dumping the afterload model data
     */
    std::string d_data_directory_name;

    /*!
     * \brief Constructor
     */
    CirculationModel(const string& object_name,
                     Pointer<Database> input_db,
                     const std::vector<Afterload_Parms> afterload_parms,
                     bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel();

    /*!
     * \brief functions for updating Windkessel variables
     */
    void windkessel_be_update(double& P, const int& src_id, const double& Q, const double& t, const double& dt);

    /*!
     * \brief Advance time-dependent data.
     */
    void advanceTimeDependentData(const double dt,
                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int U_idx,
                                  const int P_idx,
                                  const int wgt_cc_idx,
                                  const int wgt_sc_idx);

    /*!
     * \name Implementation of Serializable interface.
     */

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database point must be non-null.
     */
    void putToDatabase(Pointer<Database> db);

private:
    /*!
     * \brief Windkessel parameters
     */
    vector<double> d_R_P;
    vector<double> d_R_C;
    vector<double> d_C;
    vector<double> d_Q_sink;
    vector<double> d_P_ramp;
    vector<double> d_t_ramp;

    /*!
     * prescribed flow function for testing circ model implementation
     */
    double flow_function(double time);
    
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel(const CirculationModel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel& operator=(const CirculationModel& that);

    /*!
     * Write out source/sink state data to disk.
     */
    void writeDataFile() const;

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    void getFromRestart();
};

#endif // four_chambered_heart_circulation_model_h
