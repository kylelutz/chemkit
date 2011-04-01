/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include "forcefield.h"

#include "atom.h"
#include "foreach.h"
#include "molecule.h"
#include "constants.h"
#include "pluginmanager.h"
#include "forcefieldatom.h"
#include "forcefieldcalculation.h"

namespace chemkit {

namespace {

Float mapEnergy(const ForceFieldCalculation *calculation)
{
    return calculation->energy();
}

void reduceEnergy(Float &result, const Float &energy)
{
    result += energy;
}

} // end anonymous namespace

// === ForceFieldPrivate =================================================== //
class ForceFieldPrivate
{
    public:
        std::string name;
        ForceField::Flags flags;
        QList<ForceFieldAtom *> atoms;
        QList<ForceFieldCalculation *> calculations;
        QList<const Molecule *> molecules;
        std::string parameterSet;
        std::string parameterFile;
        std::map<std::string, std::string> parameterSets;
        std::string errorString;
};

// === ForceField ========================================================== //
/// \class ForceField forcefield.h chemkit/forcefield.h
/// \ingroup chemkit
/// \brief The ForceField class provides a generic interface to
///        molecular mechanics force fields.
///
/// The following force fields are supported in chemkit:
///     - \c amber
///     - \c mmff
///     - \c opls
///     - \c uff
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
}

/// Destroys a force field.
ForceField::~ForceField()
{
    // delete all calculations
    Q_FOREACH(ForceFieldCalculation *calculation, d->calculations){
        delete calculation;
    }

    // delete all atoms
    Q_FOREACH(ForceFieldAtom *atom, d->atoms){
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
void ForceField::setFlags(Flags flags)
{
    d->flags = flags;
}

/// Returns the flags for the force field.
ForceField::Flags ForceField::flags() const
{
    return d->flags;
}

/// Returns the number of atoms in the force field.
int ForceField::size() const
{
    return atomCount();
}

/// Returns a list of all the atoms in the force field.
QList<ForceFieldAtom *> ForceField::atoms() const
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
    return d->atoms.value(index, 0);
}

/// Returns the force field atom that represents atom.
ForceFieldAtom* ForceField::atom(const Atom *atom) const
{
    Q_FOREACH(ForceFieldAtom *forceFieldAtom, d->atoms){
        if(forceFieldAtom->atom() == atom){
            return forceFieldAtom;
        }
    }

    return 0;
}

// --- Setup --------------------------------------------------------------- //
/// Adds a molecule to the force field.
void ForceField::addMolecule(const Molecule *molecule)
{
    d->molecules.append(molecule);
}

/// Removes a molecule from the force field.
void ForceField::removeMolecule(const Molecule *molecule)
{
    d->molecules.removeAll(molecule);
}

/// Returns a list of all the molecules in the force field.
QList<const Molecule *> ForceField::molecules() const
{
    return d->molecules;
}

/// Returns the number of molecules in the force field.
int ForceField::moleculeCount() const
{
    return d->molecules.size();
}

void ForceField::addAtom(ForceFieldAtom *atom)
{
    d->atoms.append(atom);
}

void ForceField::removeAtom(ForceFieldAtom *atom)
{
    d->atoms.removeOne(atom);
}

/// Removes all of the molecules in the force field.
void ForceField::clear()
{
    Q_FOREACH(const Molecule *molecule, d->molecules){
        removeMolecule(molecule);
    }

    Q_FOREACH(ForceFieldCalculation *calculation, d->calculations){
        removeCalculation(calculation);
    }
}

/// Sets up the force field. Returns false if the setup failed.
bool ForceField::setup()
{
    return false;
}

/// Returns \c true if the force field is setup.
bool ForceField::isSetup() const
{
    Q_FOREACH(const ForceFieldCalculation *calculation, d->calculations){
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
    d->calculations.append(calculation);
}

void ForceField::removeCalculation(ForceFieldCalculation *calculation)
{
    d->calculations.removeOne(calculation);
    delete calculation;
}

/// Returns a list of all the calculations in the force field.
QList<ForceFieldCalculation *> ForceField::calculations() const
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
Float ForceField::energy() const
{
    const int parallelThreshold = 5000;

    Float energy = 0;

    if(d->calculations.size() < parallelThreshold){
        // calculate energy sequentially
        Q_FOREACH(const ForceFieldCalculation *calculation, d->calculations){
            energy += calculation->energy();
        }
    }
    else{
        // calculate energy in parallel
        energy = QtConcurrent::blockingMappedReduced(d->calculations, mapEnergy, reduceEnergy);
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
    if(d->flags.testFlag(AnalyticalGradient)){
        std::vector<Vector3> gradient(atomCount());

        Q_FOREACH(const ForceFieldCalculation *calculation, d->calculations){
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
        Float eI = atom->energy();
        Float epsilon = 1.0e-10;

        atom->moveBy(epsilon, 0, 0);
        Float eF_x = atom->energy();

        atom->moveBy(-epsilon, epsilon, 0);
        Float eF_y = atom->energy();

        atom->moveBy(0, -epsilon, epsilon);
        Float eF_z = atom->energy();

        // restore initial position
        atom->moveBy(0, 0, -epsilon);

        Float dx = (eF_x - eI) / epsilon;
        Float dy = (eF_y - eI) / epsilon;
        Float dz = (eF_z - eI) / epsilon;

        gradient[i] = Vector3(dx, dy, dz);
    }

    return gradient;
}

/// Returns the magnitude of the largest gradient.
Float ForceField::largestGradient() const
{
    if(!size()){
        return 0;
    }

    Float largest = 0;

    std::vector<Vector3> gradient = this->gradient();

    for(unsigned int i = 0; i < gradient.size(); i++){
        Float length = gradient[i].length();

        if(length > largest)
            largest = length;
    }

    return largest;
}

/// Returns the root mean square gradient.
Float ForceField::rootMeanSquareGradient() const
{
    if(!size()){
        return 0;
    }

    Float sum = 0;

    std::vector<Vector3> gradient = this->gradient();

    for(unsigned int i = 0; i < gradient.size(); i++){
        sum += gradient[i].lengthSquared();
    }

    return sqrt(sum / (3.0 * size()));
}

// --- Coordinates --------------------------------------------------------- //
/// Updates the coordinates of molecule in the force field.
void ForceField::readCoordinates(const Molecule *molecule)
{
    Q_FOREACH(const Atom *atom, molecule->atoms()){
        readCoordinates(atom);
    }
}

/// Updates the coordinates of atom in the force field.
void ForceField::readCoordinates(const Atom *atom)
{
    ForceFieldAtom *forceFieldAtom = this->atom(atom);

    if(forceFieldAtom){
        forceFieldAtom->setPosition(atom->position());
    }
}

/// Writes the coordinates to molecule from the force field.
void ForceField::writeCoordinates(Molecule *molecule) const
{
    Q_FOREACH(Atom *atom, molecule->atoms()){
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
bool ForceField::minimizationStep(Float converganceValue)
{
    // calculate gradient
    std::vector<Vector3> gradient = this->gradient();

    // perform line search
    std::vector<Point3> initialPositions(atomCount());

    Float step = 0.05;
    Float stepConv = 1e-5;
    int stepCount = 10;

    Float initialEnergy = energy();

    for(int i = 0; i < stepCount; i++){
        for(int atomIndex = 0; atomIndex < atomCount(); atomIndex++){
            ForceFieldAtom *atom = d->atoms[atomIndex];

            initialPositions[atomIndex] = atom->position();
            atom->moveBy(-gradient[atomIndex] * step);
        }

        Float finalEnergy = energy();

        // if the final energy is NaN then most likely the
        // simulation exploded so we reset the initial atom
        // positions and then 'wiggle' each atom by one
        // Angstrom in a random direction
        if(qIsNaN(finalEnergy)){
            for(int atomIndex = 0; atomIndex < atomCount(); atomIndex++){
                d->atoms[atomIndex]->setPosition(initialPositions[atomIndex]);
                d->atoms[atomIndex]->moveBy(Vector3::randomUnitVector());
            }

            // recalculate gradient
            gradient = this->gradient();

            // continue to next step
            continue;
        }

        if(finalEnergy < initialEnergy && qAbs(finalEnergy - initialEnergy) < stepConv){
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
    return rootMeanSquareGradient() < converganceValue;
}

QFuture<bool> ForceField::minimizationStepAsync(Float converganceValue)
{
    return QtConcurrent::run(this, &ForceField::minimizationStep, converganceValue);
}

// --- Geometry ------------------------------------------------------------ //
Float ForceField::distance(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return Point3::distance(a->position(), b->position());
}

Float ForceField::bondAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return bondAngleRadians(a, b, c) * chemkit::constants::RadiansToDegrees;
}

Float ForceField::bondAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return Point3::angleRadians(a->position(), b->position(), c->position());
}

Float ForceField::torsionAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return torsionAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

Float ForceField::torsionAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point3::torsionAngleRadians(a->position(), b->position(), c->position(), d->position());
}

Float ForceField::wilsonAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return wilsonAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

Float ForceField::wilsonAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point3::wilsonAngleRadians(a->position(), b->position(), c->position(), d->position());
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string that describes the last error that occured.
void ForceField::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
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
