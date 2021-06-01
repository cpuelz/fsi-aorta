// Filename: VelocityBcCoefs.h
// Created on 04 May 2007 by Boyce Griffith

#ifndef four_chambered_heart_velocity_bc_coefs_h
#define four_chambered_heart_velocity_bc_coefs_h

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// PETSC INCLUDES
#include <petscsys.h>

// VALVE TESTER INCLUDES
#include "CirculationModel.h"

// SAMRAI INCLUDES
#include <RobinBcCoefStrategy.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs : public RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs(const CirculationModel* circ_model, const int comp_idx);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                    Pointer<ArrayData<NDIM, double> >& bcoef_data,
                    Pointer<ArrayData<NDIM, double> >& gcoef_data,
                    const Pointer<hier::Variable<NDIM> >& variable,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
    IntVector<NDIM> numberOfExtensionsFillable() const;

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VelocityBcCoefs(const VelocityBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs& operator=(const VelocityBcCoefs& that);

    const CirculationModel* const d_circ_model;
    const int d_comp_idx;
};

#endif // four_chambered_heart_velocity_bc_coefs_h
