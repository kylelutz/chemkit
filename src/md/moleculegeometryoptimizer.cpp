/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "moleculegeometryoptimizer.h"

#include <boost/make_shared.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/concurrent.h>
#include <chemkit/cartesiancoordinates.h>

#include "forcefield.h"
#include "integrator.h"

namespace chemkit {

namespace {

class SteepestDescentIntegrator : public Integrator
{
public:
    void integrate() CHEMKIT_OVERRIDE;
};

void SteepestDescentIntegrator::integrate()
{
    CartesianCoordinates *coordinates = this->coordinates();
    boost::shared_ptr<Potential> potential = this->potential();

    if(!potential || !coordinates){
        return;
    }

    // optimization parameters
    Real step = 0.05;
    Real stepConv = 1e-5;
    size_t stepCount = 10;

    // calculate initial energy and gradient
    Real initialEnergy = potential->energy(coordinates);
    std::vector<Vector3> gradient = potential->gradient(coordinates);

    // perform line search
    for(size_t i = 0; i < stepCount; i++){
        // save initial coordinates
        CartesianCoordinates initialCoordinates = *coordinates;

        // move each atom against its gradient
        for(size_t atomIndex = 0; atomIndex < potential->size(); atomIndex++){
            (*coordinates)[atomIndex] += -gradient[atomIndex] * step;
        }

        // calculate new energy
        Real finalEnergy = potential->energy(coordinates);

        // if the final energy is NaN then most likely the
        // simulation exploded so we reset the initial atom
        // positions and then 'wiggle' each atom by one
        // Angstrom in a random direction
        if((boost::math::isnan)(finalEnergy)){
            for(size_t atomIndex = 0; atomIndex < potential->size(); atomIndex++){
                Point3 position = initialCoordinates.position(atomIndex);
                position += Vector3::Random().normalized();
                coordinates->setPosition(atomIndex, position);
            }

            // recalculate gradient
            gradient = potential->gradient(coordinates);

            // continue to next step
            continue;
        }

        if(finalEnergy < initialEnergy && std::abs(finalEnergy - initialEnergy) < stepConv){
            break;
        }
        else if(finalEnergy < initialEnergy){
            // we reduced the energy, so set a bigger step size
            step *= 2;

            // maximum step size is 1
            if(step > 1){
                step = 1;
            }

            // the initial energy for the next step
            // is the final energy of this step
            initialEnergy = finalEnergy;
        }
        else if(finalEnergy > initialEnergy){
            // we went too far, so reset initial atom positions
            *coordinates = initialCoordinates;

            // and reduce step size
            step *= 0.1;
        }
    }
}

} // end anonymous namespace

// === MoleculeGeometryOptimizerPrivate ==================================== //
class MoleculeGeometryOptimizerPrivate
{
public:
    Molecule *molecule;
    boost::shared_ptr<ForceField> forceField;
    std::string forceFieldName;
    std::string errorString;
    boost::shared_ptr<Integrator> integrator;
};

// === MoleculeGeometryOptimizer =========================================== //
/// \class MoleculeGeometryOptimizer moleculegeometryoptimizer.h chemkit/moleculegeometryoptimizer.h
/// \ingroup chemkit-md
/// \brief The MoleculeGeometryOptimizer class performs geometry
///        optimization for a single molecule.
///
/// The MoleculeGeometryOptimizer class is a convenience class meant
/// to simplify the process of setting up a force field and
/// performing an energy minimization run for a single molecule.
///
/// By default the UFF force field is used.
///
/// The easiest way to optimize the geometry for a molecule is to
/// use the optimizeCoordinate() static method as follows:
/// \code
/// MoleculeGeometryOptimizer::optimizeCoordinates(molecule);
/// \endcode
///
/// This class along with the CoordinatePredictor can be used to
/// generate 3D coordinates for a molecule. The following example
/// shows how to create a phenol molecule from its SMILES formula
/// and generate a rough set of 3D coordinates.
/// \code
/// // create phenol molecule from its formula
/// Molecule phenol("c1ccccc1O", "smiles");
///
/// // predict an initial set of 3D coordinate
/// CoordinatePredictor::predictCoordinates(&phenol);
///
/// // optimize the predicted coordinates
/// MoleculeGeometryOptimizer::optimizeCoordinates(&phenol);
/// \endcode
///
/// \see ForceField

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new geometry optimizer for \p molecule.
MoleculeGeometryOptimizer::MoleculeGeometryOptimizer(Molecule *molecule)
    : d(new MoleculeGeometryOptimizerPrivate)
{
    d->molecule = molecule;
    d->forceFieldName = "uff";
    d->integrator = boost::make_shared<SteepestDescentIntegrator>();
}

/// Destroys the geometry optmizer object.
MoleculeGeometryOptimizer::~MoleculeGeometryOptimizer()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the geometry optimizer to \p molecule.
void MoleculeGeometryOptimizer::setMolecule(Molecule *molecule)
{
    if(molecule != d->molecule){
        d->molecule = molecule;
    }
}

/// Returns the molecule for the geometry optimizer.
Molecule* MoleculeGeometryOptimizer::molecule() const
{
    return d->molecule;
}

/// Sets the force field to use to \p forceField. Returns \c false
/// if \p forceField is not supported.
///
/// \see ForceField
bool MoleculeGeometryOptimizer::setForceField(const std::string &forceField)
{
    d->forceFieldName = forceField;

    return true;
}

/// Returns the name of the force field being used for geometry
/// optimization.
std::string MoleculeGeometryOptimizer::forceField() const
{
    return d->forceFieldName;
}

// --- Energy -------------------------------------------------------------- //
/// Returns the current energy of the force field.
Real MoleculeGeometryOptimizer::energy() const
{
    if(!d->forceField){
        return 0;
    }

    return d->integrator->energy();
}

// --- Optimization -------------------------------------------------------- //
/// Sets up the force field. Returns \c false if an error occurred.
bool MoleculeGeometryOptimizer::setup()
{
    if(!d->molecule){
        d->errorString = "No molecule specified";
        return false;
    }

    d->forceField = boost::shared_ptr<ForceField>(ForceField::create(d->forceFieldName));
    if(!d->forceField){
        d->errorString = "Force field '" + d->forceFieldName + "' is not supported.";
        return false;
    }

    d->integrator->setPotential(d->forceField);
    d->integrator->setCoordinates(d->molecule->coordinates());

    d->forceField->setTopologyFromMolecule(d->molecule);
    if(!d->forceField->setup()){
        d->errorString = "Failed to setup force field.";
        return false;
    }

    return true;
}

/// Performs a single geometry optimization step. Returns \c true if
/// converged. The minimization is considered converged when the
/// root mean square gradient is below the convergance value.
void MoleculeGeometryOptimizer::step()
{
    if(!d->molecule || !d->forceField){
        return;
    }

    // perform a single integration step
    d->integrator->integrate();
}

/// Returns \c true if the optimization algorithm has converged. By
/// default, the algorithm is considered converged when the
/// root-mean-square gradient of the force field falls below \c 0.1.
bool MoleculeGeometryOptimizer::converged()
{
    if(!d->forceField){
        return false;
    }

    // check for convergance
    return d->integrator->rmsg() < 0.1;
}

/// Optimizes the geometry of the molecule. Returns \c true if the
/// optimization algorithm converged.
bool MoleculeGeometryOptimizer::optimize()
{
    if(!setup()){
        return false;
    }

    while(!converged()){
        step();
    }

    // write the optimized coordinates to the molecule
    writeCoordinates();

    return true;
}

/// Writes the optimized coordinates to the molecule.
void MoleculeGeometryOptimizer::writeCoordinates()
{
    if(!d->molecule){
        return;
    }

    const CartesianCoordinates *coordinates = d->integrator->coordinates();

    for(size_t i = 0; i < d->molecule->size(); i++){
        d->molecule->atom(i)->setPosition(coordinates->position(i));
    }
}

// --- Error Handling ------------------------------------------------------ //
/// Returns a string describing the last error that occurred.
std::string MoleculeGeometryOptimizer::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Optimizes the geometry of \p molecule.
bool MoleculeGeometryOptimizer::optimizeCoordinates(Molecule *molecule)
{
    return MoleculeGeometryOptimizer(molecule).optimize();
}

/// Runs the optimizeCoordinates() method asynchronously and returns
/// a future containing the result.
///
/// \internal
boost::shared_future<bool> MoleculeGeometryOptimizer::optimizeCoordinatesAsync(Molecule *molecule)
{
  return chemkit::concurrent::run(
      boost::bind(&MoleculeGeometryOptimizer::optimizeCoordinates, molecule));
}

} // end chemkit namespace
