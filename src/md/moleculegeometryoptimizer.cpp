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

#include <boost/math/special_functions/fpclassify.hpp>

#include "forcefield.h"

namespace chemkit {

// === MoleculeGeometryOptimizerPrivate ==================================== //
class MoleculeGeometryOptimizerPrivate
{
public:
    Molecule *molecule;
    ForceField *forceField;
    std::string forceFieldName;
    std::string errorString;
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
    d->forceField = 0;
    d->forceFieldName = "uff";
}

/// Destroys the geometry optmizer object.
MoleculeGeometryOptimizer::~MoleculeGeometryOptimizer()
{
    delete d->forceField;
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

    return d->forceField->energy();
}

// --- Optimization -------------------------------------------------------- //
/// Sets up the force field. Returns \c false if an error occurred.
bool MoleculeGeometryOptimizer::setup()
{
    if(!d->molecule){
        d->errorString = "No molecule specified";
        return false;
    }

    delete d->forceField;

    d->forceField = chemkit::ForceField::create(d->forceFieldName);
    if(!d->forceField){
        d->errorString = "Force field '" + d->forceFieldName + "' is not supported.";
        return false;
    }

    d->forceField->setMolecule(d->molecule);

    if(!d->forceField->setup()){
        d->errorString = "Failed to setup force field.";
        return false;
    }

    return true;
}

/// Performs a single geometry optimization step. Returns \c true if
/// converged. The minimization is considered converged when the
/// root mean square gradient is below the convergance value.
bool MoleculeGeometryOptimizer::step()
{
    if(!d->molecule || !d->forceField){
        return false;
    }

    Real converganceValue = 0.1;

    // calculate gradient
    std::vector<Vector3> gradient = d->forceField->gradient();

    // perform line search
    std::vector<Point3> initialPositions(d->forceField->atomCount());

    Real step = 0.05;
    Real stepConv = 1e-5;
    int stepCount = 10;

    Real initialEnergy = d->forceField->energy();

    for(int i = 0; i < stepCount; i++){
        for(int atomIndex = 0; atomIndex < d->forceField->atomCount(); atomIndex++){
            ForceFieldAtom *atom = d->forceField->atom(atomIndex);

            initialPositions[atomIndex] = atom->position();
            atom->moveBy(-gradient[atomIndex] * step);
        }

        Real finalEnergy = d->forceField->energy();

        // if the final energy is NaN then most likely the
        // simulation exploded so we reset the initial atom
        // positions and then 'wiggle' each atom by one
        // Angstrom in a random direction
        if((boost::math::isnan)(finalEnergy)){
            for(int atomIndex = 0; atomIndex < d->forceField->atomCount(); atomIndex++){
                ForceFieldAtom *atom = d->forceField->atom(atomIndex);

                atom->setPosition(initialPositions[atomIndex]);
                atom->moveBy(Vector3::Random().normalized());
            }

            // recalculate gradient
            gradient = d->forceField->gradient();

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
            for(int atomIndex = 0; atomIndex < d->forceField->atomCount(); atomIndex++){
                ForceFieldAtom *atom = d->forceField->atom(atomIndex);

                atom->setPosition(initialPositions[atomIndex]);
            }

            // and reduce step size
            step *= 0.1;
        }
    }

    // check for convergance
    return d->forceField->rmsg() < converganceValue;
}

/// Optimizes the geometry of the molecule. Returns \c true if the
/// optimization algorithm converged.
bool MoleculeGeometryOptimizer::optimize()
{
    if(!setup()){
        return false;
    }

    bool done = false;
    while(!done){
        done = step();
    }

    d->forceField->writeCoordinates(d->molecule);

    return done;
}

/// Writes the optimized coordinates to the molecule.
void MoleculeGeometryOptimizer::writeCoordinates()
{
    if(!d->molecule || !d->forceField){
        return;
    }

    d->forceField->writeCoordinates(d->molecule);
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

} // end chemkit namespace
