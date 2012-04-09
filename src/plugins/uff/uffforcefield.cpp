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

// The uff plugin implements the UFF force field.
//
// See:
//   - [Rappe 1992]
//   - http://towhee.sourceforge.net/forcefields/uff.html

#include "uffforcefield.h"

#include <boost/algorithm/string.hpp>

#include "uffatomtyper.h"
#include "uffparameters.h"
#include "uffcalculation.h"

#include <chemkit/foreach.h>
#include <chemkit/topology.h>

// --- Construction and Destruction ---------------------------------------- //
UffForceField::UffForceField()
    : chemkit::ForceField("uff")
{
    m_parameters = new UffParameters;

    setFlags(chemkit::ForceField::AnalyticalGradient);
}

UffForceField::~UffForceField()
{
    delete m_parameters;
}

// --- Parameters ---------------------------------------------------------- //
const UffParameters* UffForceField::parameters() const
{
    return m_parameters;
}

// --- Setup --------------------------------------------------------------- //
bool UffForceField::setup()
{
    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();
    if(!topology){
        return false;
    }

    // bond strech
    foreach(const chemkit::Topology::BondedInteraction &interaction, topology->bondedInteractions()){
        addCalculation(new UffBondStrechCalculation(interaction[0],
                                                    interaction[1]));
    }

    // angle bend
    foreach(const chemkit::Topology::AngleInteraction &interaction, topology->angleInteractions()){
        addCalculation(new UffAngleBendCalculation(interaction[0],
                                                   interaction[1],
                                                   interaction[2]));
    }

    // torsion
    foreach(const chemkit::Topology::TorsionInteraction &interaction, topology->torsionInteractions()){
        addCalculation(new UffTorsionCalculation(interaction[0],
                                                 interaction[1],
                                                 interaction[2],
                                                 interaction[3]));
    }

    // inversion
    foreach(const chemkit::Topology::ImproperTorsionInteraction &interaction, topology->improperTorsionInteractions()){
        // type for the center atom
        std::string typeB = topology->type(interaction[1]);

        if(boost::starts_with(typeB, "C_") ||
           boost::starts_with(typeB, "N_") ||
           boost::starts_with(typeB, "P_") ||
           boost::starts_with(typeB, "As") ||
           boost::starts_with(typeB, "Sb") ||
           boost::starts_with(typeB, "Bi")){
            addCalculation(new UffInversionCalculation(interaction[0],
                                                       interaction[1],
                                                       interaction[2],
                                                       interaction[3]));
            addCalculation(new UffInversionCalculation(interaction[0],
                                                       interaction[1],
                                                       interaction[3],
                                                       interaction[2]));
            addCalculation(new UffInversionCalculation(interaction[2],
                                                       interaction[1],
                                                       interaction[0],
                                                       interaction[3]));
        }
    }

    // van der waals
    foreach(const chemkit::Topology::NonbondedInteraction &interaction, topology->nonbondedInteractions()){
        addCalculation(new UffVanDerWaalsCalculation(interaction[0],
                                                     interaction[1]));
    }

    bool ok = true;

    foreach(chemkit::ForceFieldCalculation *calculation, calculations()){
        bool setup = static_cast<UffCalculation *>(calculation)->setup();

        if(!setup){
            ok = false;
        }

        setCalculationSetup(calculation, setup);
    }

    return ok;
}

/// Returns \c true if \p atom is in group six of the periodic table.
bool UffForceField::isGroupSix(size_t atom) const
{
    std::string type = topology()->type(atom);

    return boost::starts_with(type, "O_") ||
           boost::starts_with(type, "S_") ||
           boost::starts_with(type, "Se") ||
           boost::starts_with(type, "Te") ||
           boost::starts_with(type, "Po");
}
