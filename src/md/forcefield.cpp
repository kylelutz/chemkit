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

#include <chemkit/foreach.h>
#include <chemkit/constants.h>
#include <chemkit/concurrent.h>
#include <chemkit/pluginmanager.h>
#include <chemkit/cartesiancoordinates.h>

#include "topology.h"
#include "topologybuilder.h"
#include "forcefieldcalculation.h"

namespace chemkit {

// === ForceFieldPrivate =================================================== //
class ForceFieldPrivate
{
public:
    std::string name;
    int flags;
    boost::shared_ptr<Topology> topology;
    std::vector<ForceFieldCalculation *> calculations;
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
/// molecule using the UFF force field.
///
/// \code
/// // create the uff force field
/// ForceField *forceField = ForceField::create("uff");
///
/// // set the topology for the force field
/// forceField->setTopologyFromMolecule(molecule);
///
/// // setup the force field
/// forceField->setup();
///
/// // calculate the total energy
/// double energy = forceField->energy(molecule->coordinates());
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
ForceField::ForceField(const std::string &name)
    : d(new ForceFieldPrivate)
{
    d->name = name;
    d->flags = 0;
}

/// Destroys a force field.
ForceField::~ForceField()
{
    // delete all calculations
    foreach(ForceFieldCalculation *calculation, d->calculations){
        delete calculation;
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
size_t ForceField::size() const
{
    if(!d->topology){
        return 0;
    }

    return d->topology->size();
}

// --- Setup --------------------------------------------------------------- //
/// Sets the topology for the force field to \p topology.
void ForceField::setTopology(const boost::shared_ptr<Topology> &topology)
{
    d->topology = topology;

    // remove old calculations
    foreach(ForceFieldCalculation *calculation, d->calculations){
        delete calculation;
    }
    d->calculations.clear();
}

/// Builds a topology for the molecule and sets it with setTopology().
///
/// \see TopologyBuilder
void ForceField::setTopologyFromMolecule(const Molecule *molecule)
{
    TopologyBuilder builder;
    builder.setAtomTyper(name());
    builder.setPartialChargeModel(name());
    builder.addMolecule(molecule);
    setTopology(builder.topology());
}

/// Returns the topology for the force field.
boost::shared_ptr<Topology> ForceField::topology() const
{
    return d->topology;
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
    calculation->setForceField(this);

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
size_t ForceField::calculationCount() const
{
    return d->calculations.size();
}

void ForceField::setCalculationSetup(ForceFieldCalculation *calculation, bool setup)
{
    calculation->setSetup(setup);
}

/// \copydoc Potential::energy()
Real ForceField::energy(const CartesianCoordinates *coordinates) const
{
    Real energy = 0;

    foreach(const ForceFieldCalculation *calculation, d->calculations){
        energy += calculation->energy(coordinates);
    }

    return energy;
}

/// \copydoc Potential::gradient()
std::vector<Vector3> ForceField::gradient(const CartesianCoordinates *coordinates) const
{
    if(d->flags & AnalyticalGradient){
        std::vector<Vector3> gradient(size());
        std::fill(gradient.begin(), gradient.end(), Vector3(0, 0, 0));

        foreach(const ForceFieldCalculation *calculation, d->calculations){
            std::vector<Vector3> atomGradients = calculation->gradient(coordinates);

            for(size_t i = 0; i < atomGradients.size(); i++){
                gradient[calculation->atom(i)] += atomGradients[i];
            }
        }

        return gradient;
    }
    else{
        return numericalGradient(coordinates);
    }
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
