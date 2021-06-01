// Filename: SourceDistributer.h
// Created on 02 Jun 2018 by Boyce Griffith

#ifndef four_chambered_heart_source_distributer_h
#define four_chambered_heart_source_distributer_h

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// HEART MODEL INCLUDES
#include "Source.h"

// IBAMR INCLUDES
#include <ibamr/INSHierarchyIntegrator.h>

// STDLIB INCLUDES
#include <vector>

// NAMESPACE
#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class SourceDistributer is an implementation of the strategy class
 * CartGridFunction that is used to specify source/sink distributions.
 */
class SourceDistributer : public CartGridFunction
{
public:
    /*!
     * \brief Constructor
     */
    SourceDistributer(std::vector<std::unique_ptr<Source> > sources,
                      const INSHierarchyIntegrator* fluid_solver,
                      Pointer<PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Destructor.
     */
    virtual ~SourceDistributer();

    /*!
     * \name Implementation of CartGridFunction interface.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Set data on the specified patch interior.
     */
    void setDataOnPatch(int data_idx,
                        Pointer<hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    /*!
     * \brief Interpolate the pressure to the sources.
     */
    void interpolatePressure(int p_data_idx, double data_time);

    /*!
     * \brief Update the sources with Eulerian data.
     */
    void updateSources(IBFEMethod* ib_method_ops, const double loop_time, const double dt);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SourceDistributer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SourceDistributer(const SourceDistributer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SourceDistributer& operator=(const SourceDistributer& that);

    std::vector<std::unique_ptr<Source> > d_sources;
    const INSHierarchyIntegrator* const d_fluid_solver;
    Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <SourceDistributer.I>

//////////////////////////////////////////////////////////////////////////////

#endif // four_chambered_heart_source_distributer_h
