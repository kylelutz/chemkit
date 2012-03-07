/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "forcefield.h"

#include <boost/math/special_functions/fpclassify.hpp>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/geometry.h>
#include <chemkit/molecule.h>
#include <chemkit/constants.h>
#include <chemkit/pluginmanager.h>

#include "forcefieldatom.h"
#include "forcefieldcalculation.h"

namespace chemkit {

// === ForceFieldPrivate =================================================== //
class ForceFieldPrivate
{
public:
    std::string name;
    int flags;
    std::vector<ForceFieldAtom *> atoms;
    std::vector<ForceFieldCalculation *> calculations;
    const Molecule *molecule;
    std::string parameterSet;
    std::string parameterFile;
    std::map<std::string, std::string> parameterSets;
    std::string errorString;
};

// === ForceField ========================================================== //
/// \class ForceField forcefield.h chemkit/forcefield.h
/// \ingroup chemkit-md
/// \brief The ForceField class provides a generic interface to
///        molecular mechanics force fields.
///
/// A list of supported force fields is available at:
/// http://wiki.chemkit.org/Features#Force_Fields
///
/// The following example shows how to calculate the energy of a
/// molecule using the uff force field.
///
/// \code
/// // create the uff force field
/// ForceField *forceField = ForceField::create("uff");
///
/// // add the molecule to the force field
/// forceField->addMolecule(molecule);
///
/// // setup the force field
/// forceField->setup();
///
/// // calculate the total energy
/// Float energy = forceField->energy();
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
ForceField::ForceField(const std::string &name)
    : d(new ForceFieldPrivate)
{
    d->name = name;
    d->molecule = 0;
    d->flags = 0;
}

/// Destroys a force field.
ForceField::~ForceField()
{
    // delete all calculations
    foreach(ForceFieldCalculation *calculation, d->calculations){
        delete calculation;
    }

    // delete all atoms
    foreach(ForceFieldAtom *atom, d->atoms){
        delete atom;
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the force field.
std::string ForceField::name() const
{
    return d->name;
}

/// Sets the flags for the force field to \p flags.
void ForceField::setFlags(int flags)
{
    d->flags = flags;
}

/// Returns the flags for the force field.
int ForceField::flags() const
{
    return d->flags;
}

/// Returns the number of atoms in the force field.
int ForceField::size() const
{
    return atomCount();
}

/// Returns a list of all the atoms in the force field.
std::vector<ForceFieldAtom *> ForceField::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the force field.
int ForceField::atomCount() const
{
    return d->atoms.size();
}

/// Returns the atom at index.
ForceFieldAtom* ForceField::atom(int index) const
{
    return d->atoms[index];
}

/// Returns the force field atom that represents atom.
ForceFieldAtom* ForceField::atom(const Atom *atom) const
{
    foreach(ForceFieldAtom *forceFieldAtom, d->atoms){
        if(forceFieldAtom->atom() == atom){
            return forceFieldAtom;
        }
    }

    return 0;
}

// --- Setup --------------------------------------------------------------- //
/// Sets the molecule for the force field to \p molecule.
void ForceField::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    foreach(ForceFieldAtom *atom, d->atoms){
        delete atom;
    }
    d->atoms.clear();

    foreach(ForceFieldCalculation *calculation, d->calculations){
        delete calculation;
    }
    d->calculations.clear();
}

/// Returns the molecule for the force field.
const Molecule* ForceField::molecule() const
{
    return d->molecule;
}

void ForceField::addAtom(ForceFieldAtom *atom)
{
    d->atoms.push_back(atom);
}

void ForceField::removeAtom(ForceFieldAtom *atom)
{
    d->atoms.erase(std::remove(d->atoms.begin(), d->atoms.end(), atom));
}

/// Removes all of the molecules in the force field.
void ForceField::clear()
{
    setMolecule(0);
}

/// Sets up the force field. Returns false if the setup failed.
bool ForceField::setup()
{
    return false;
}

/// Returns \c true if the force field is setup.
bool ForceField::isSetup() const
{
    foreach(const ForceFieldCalculation *calculation, d->calculations){
        if(!calculation->isSetup()){
            return false;
        }
    }

    return true;
}

// --- Parameters ---------------------------------------------------------- //
void ForceField::addParameterSet(const std::string &name, const std::string &fileName)
{
    d->parameterSets[name] = fileName;
}

void ForceField::removeParameterSet(const std::string &name)
{
    d->parameterSets.erase(name);
}

void ForceField::setParameterSet(const std::string &name)
{
    std::map<std::string, std::string>::iterator element = d->parameterSets.find(name);
    if(element == d->parameterSets.end()){
        return;
    }

    d->parameterSet = name;
    d->parameterFile = element->second;
}

std::string ForceField::parameterSet() const
{
    return d->parameterSet;
}

std::vector<std::string> ForceField::parameterSets() const
{
    std::vector<std::string> parameterSets;

    std::pair<std::string, std::string> element;
    foreach(element, d->parameterSets){
        parameterSets.push_back(element.first);
    }

    return parameterSets;
}

void ForceField::setParameterFile(const std::string &fileName)
{
    d->parameterFile = fileName;
}

std::string ForceField::parameterFile() const
{
    return d->parameterFile;
}

// --- Calculations -------------------------------------------------------- //
void ForceField::addCalculation(ForceFieldCalculation *calculation)
{
    d->calculations.push_back(calculation);
}

void ForceField::removeCalculation(ForceFieldCalculation *calculation)
{
    d->calculations.erase(std::remove(d->calculations.begin(), d->calculations.end(), calculation));
    delete calculation;
}

/// Returns a list of all the calculations in the force field.
std::vector<ForceFieldCalculation *> ForceField::calculations() const
{
    return d->calculations;
}

/// Returns the number of calculations in the force field.
int ForceField::calculationCount() const
{
    return d->calculations.size();
}

void ForceField::setCalculationSetup(ForceFieldCalculation *calculation, bool setup)
{
    calculation->setSetup(setup);
}

/// Calculates and returns the total energy of the system. Energy is
/// in kcal/mol. If the force field is not setup this method will
/// return \c 0.
Real ForceField::energy() const
{
    Real energy = 0;

    foreach(const ForceFieldCalculation *calculation, d->calculations){
        energy += calculation->energy();
    }

    return energy;
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom in the force field.
///
/** \f[ \nabla E = \left[
///                \begin{array}{ccc}
///                    \frac{\partial E}{\partial x_{0}} &
///                    \frac{\partial E}{\partial y_{0}} &
///                    \frac{\partial E}{\partial z_{0}} \\
///                    \frac{\partial E}{\partial x_{1}} &
///                    \frac{\partial E}{\partial y_{1}} &
///                    \frac{\partial E}{\partial z_{1}} \\
///                    \vdots & \vdots & \vdots \\
///                    \frac{\partial E}{\partial x_{n}} &
///                    \frac{\partial E}{\partial y_{n}} &
///                    \frac{\partial E}{\partial z_{n}}
///                \end{array}
///                \right]
/// \f]
**/
std::vector<Vector3> ForceField::gradient() const
{
    if(d->flags & AnalyticalGradient){
        std::vector<Vector3> gradient(atomCount());
        std::fill(gradient.begin(), gradient.end(), Vector3(0, 0, 0));

        foreach(const ForceFieldCalculation *calculation, d->calculations){
            std::vector<Vector3> atomGradients = calculation->gradient();

            for(unsigned int i = 0; i < atomGradients.size(); i++){
                const ForceFieldAtom *atom = calculation->atom(i);

                gradient[atom->index()] += atomGradients[i];
            }
        }

        return gradient;
    }
    else{
        return numericalGradient();
    }
}

/// Returns the gradient of the energy with respect to the
/// coordinates of each atom in the force field. The gradient is
/// calculated numerically.
///
/// \see ForceField::gradient()
std::vector<Vector3> ForceField::numericalGradient() const
{
    std::vector<Vector3> gradient(atomCount());

    for(int i = 0; i < atomCount(); i++){
        ForceFieldAtom *atom = d->atoms[i];

        // initial energy
        Real eI = atom->energy();
        Real epsilon = 1.0e-10;

        atom->moveBy(epsilon, 0, 0);
        Real eF_x = atom->energy();

        atom->moveBy(-epsilon, epsilon, 0);
        Real eF_y = atom->energy();

        atom->moveBy(0, -epsilon, epsilon);
        Real eF_z = atom->energy();

        // restore initial position
        atom->moveBy(0, 0, -epsilon);

        Real dx = (eF_x - eI) / epsilon;
        Real dy = (eF_y - eI) / epsilon;
        Real dz = (eF_z - eI) / epsilon;

        gradient[i] = Vector3(dx, dy, dz);
    }

    return gradient;
}

/// Returns the root mean square gradient.
Real ForceField::rmsg() const
{
    if(!size()){
        return 0;
    }

    Real sum = 0;

    std::vector<Vector3> gradient = this->gradient();

    for(unsigned int i = 0; i < gradient.size(); i++){
        sum += gradient[i].squaredNorm();
    }

    return sqrt(sum / (3.0 * size()));
}

// --- Coordinates --------------------------------------------------------- //
/// Writes the coordinates to molecule from the force field.
void ForceField::writeCoordinates(Molecule *molecule) const
{
    foreach(Atom *atom, molecule->atoms()){
        writeCoordinates(atom);
    }
}

/// Writes the coordinates to atom from the force field.
void ForceField::writeCoordinates(Atom *atom) const
{
    const ForceFieldAtom *forceFieldAtom = this->atom(atom);

    if(forceFieldAtom){
        atom->setPosition(forceFieldAtom->position());
    }
}

// --- Energy Minimization ------------------------------------------------- //
/// Perform one step of energy minimization. Returns \c true if
/// converged. The minimization is considered converged when the
/// root mean square gradient is below \p converganceValue.
bool ForceField::minimizationStep(Real converganceValue)
{
    // calculate gradient
    std::vector<Vector3> gradient = this->gradient();

    // perform line search
    std::vector<Point3> initialPositions(atomCount());

    Real step = 0.05;
    Real stepConv = 1e-5;
    int stepCount = 10;

    Real initialEnergy = energy();

    for(int i = 0; i < stepCount; i++){
        for(int atomIndex = 0; atomIndex < atomCount(); atomIndex++){
            ForceFieldAtom *atom = d->atoms[atomIndex];

            initialPositions[atomIndex] = atom->position();
            atom->moveBy(-gradient[atomIndex] * step);
        }

        Real finalEnergy = energy();

        // if the final energy is NaN then most likely the
        // simulation exploded so we reset the initial atom
        // positions and then 'wiggle' each atom by one
        // Angstrom in a random direction
        if((boost::math::isnan)(finalEnergy)){
            for(int atomIndex = 0; atomIndex < atomCount(); atomIndex++){
                d->atoms[atomIndex]->setPosition(initialPositions[atomIndex]);
                d->atoms[atomIndex]->moveBy(Vector3::Random().normalized());
            }

            // recalculate gradient
            gradient = this->gradient();

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
            for(int atomIndex = 0; atomIndex < atomCount(); atomIndex++){
                d->atoms[atomIndex]->setPosition(initialPositions[atomIndex]);
            }

            // and reduce step size
            step *= 0.1;
        }
    }

    // check for convergance
    return rmsg() < converganceValue;
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string that describes the last error that occurred.
void ForceField::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occurred.
std::string ForceField::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Create a new force field from \p name. If \p name is invalid or
/// a force field with \p name is not available \c 0 is returned.
ForceField* ForceField::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<ForceField>(name);
}

/// Returns a list of names of all supported force fields.
std::vector<std::string> ForceField::forceFields()
{
    return PluginManager::instance()->pluginClassNames<ForceField>();
}

} // end chemkit namespace
