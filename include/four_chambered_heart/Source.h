/*
 * File:   Source.h
 * Author: cpuelz
 *
 * Created on June 2, 2018, 2:48 PM
 */

#ifndef four_chambered_heart_source_h
#define four_chambered_heart_source_h

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/Serializable.h>

// IBTK INCLUDES
#include <ibtk/ibtk_utilities.h>

#include "ibamr/IBFEMethod.h"
#include <ibamr/app_namespaces.h>

class Source : public Serializable
{
public:

    /*!
     * \brief Keep track of time for output.
     */
    double d_time;

    /*!
     * \brief file stream for output
     */
    std::ofstream d_data_stream;

    /*!
     * \brief The object name.
     */
    string d_object_name;

    /*!
     * \brief Whether the object is registered with the restart manager.
     */
    bool d_registered_for_restart;

    /*!
     * \brief whether the source is just sampling the pressure and doing nothing
     * else.
     */
    bool d_is_pressure_meter;
    
    /*!
     * \brief whether the source prescribing the flow directly
     * or if the flow is driven by the pressure gradient.
     */
    bool d_is_pressure_source;
    
     /*!
     * \brief this is the prescribed flow and is needed is d_is_pressure_source
      *  = FALSE
     */
    bool d_prescribed_flow;

    /*!
     * \brief part id.
     */
    int d_part;

    /*!
     * \brief name of sidesets from which the source centroid is computed.
     */
    SAMRAI::tbox::Array<std::string> d_sideset_names;

    /*!
     * \brief sideset ids from which the source centroid is computed.
     */
    std::vector<int> d_sideset_IDs;

    /*!
     * \brief location of source in reference configuration.
     */
    std::vector<double> d_initial_location;

    /*!
     * \brief location of source in current configurations
     */
    std::vector<double> d_current_location;

    /*!
     * \brief values of flow and pressure at the source.
     */
    double d_q_src, d_p_src;

    /*!
     * \brief where fluid source data file is stored
     */
    std::string d_data_directory_name;

    /*!
     * \brief filename containing all of the output data
     */
    std::string d_data_filename;

    // constructor
    Source(const string& object_name,
           const int part,
           Pointer<Database> source_db,
           IBFEMethod* ib_method_ops,
           bool register_for_restart = true);

    // update current location
    void updateCurrentLocation(IBFEMethod* ib_method_ops);

    // update source strength
    void updateSourceStrength(const double loop_time, const double dt, const int source_number);

    // destructor
    ~Source();

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
    // parameters for the source
    double d_R_src, d_L_src, d_P_reservoir;
    double d_ramp_time;

    // This class (like most of IBAMR) assumes that we never regrid the
    // Lagrangian data. Hence, to improve performance, store the dof indices
    // of boundary faces on the given sidesets here so that we don't need to
    // loop across the mesh to get them.
    IBAMR::IBFEMethod *d_ibfe_method = nullptr;
    std::array<std::vector<libMesh::dof_id_type>, NDIM> d_locally_owned_side_dofs;
    int d_nodes_per_side = -1;
    int d_local_num_elem_sides = -1;

    /*!
     * Write out data to disk
     */
    void writeDataFile();

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

#endif // four_chambered_heart_source_h
