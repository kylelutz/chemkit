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

#include "uffatomtyper.h"
#include "uffparameters.h"
#include "uffcalculation.h"

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/forcefieldinteractions.h>

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
    const chemkit::Molecule *molecule = this->molecule();
    if(!molecule){
        return false;
    }

    std::map<const chemkit::Atom *, chemkit::ForceFieldAtom *> atoms;

    UffAtomTyper typer(molecule);

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        chemkit::ForceFieldAtom *forceFieldAtom = new chemkit::ForceFieldAtom(this, atom);
        atoms[atom] = forceFieldAtom;
        addAtom(forceFieldAtom);
        forceFieldAtom->setType(typer.typeString(atom).c_str());
    }

    chemkit::ForceFieldInteractions interactions(molecule, this);

    // bond strech
    std::pair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> bondedPair;
    foreach(bondedPair, interactions.bondedPairs()){
        addCalculation(new UffBondStrechCalculation(bondedPair.first,
                                                    bondedPair.second));
    }

    // angle bend
    std::vector<const chemkit::ForceFieldAtom *> angleGroup;
    foreach(angleGroup, interactions.angleGroups()){
        addCalculation(new UffAngleBendCalculation(angleGroup[0],
                                                   angleGroup[1],
                                                   angleGroup[2]));
    }

    // torsion
    std::vector<const chemkit::ForceFieldAtom *> torsionGroup;
    foreach(torsionGroup, interactions.torsionGroups()){
        addCalculation(new UffTorsionCalculation(torsionGroup[0],
                                                 torsionGroup[1],
                                                 torsionGroup[2],
                                                 torsionGroup[3]));
    }

    // inversion
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->neighborCount() == 3 && (atom->is(chemkit::Atom::Carbon) ||
                                          atom->is(chemkit::Atom::Nitrogen) ||
                                          atom->is(chemkit::Atom::Phosphorus) ||
                                          atom->is(chemkit::Atom::Arsenic) ||
                                          atom->is(chemkit::Atom::Antimony) ||
                                          atom->is(chemkit::Atom::Bismuth))){
            std::vector<chemkit::Atom *> neighbors(atom->neighbors().begin(),
                                                   atom->neighbors().end());

            addCalculation(new UffInversionCalculation(atoms[neighbors[0]],
                                                       atoms[atom],
                                                       atoms[neighbors[1]],
                                                       atoms[neighbors[2]]));
            addCalculation(new UffInversionCalculation(atoms[neighbors[0]],
                                                       atoms[atom],
                                                       atoms[neighbors[2]],
                                                       atoms[neighbors[1]]));
            addCalculation(new UffInversionCalculation(atoms[neighbors[1]],
                                                       atoms[atom],
                                                       atoms[neighbors[2]],
                                                       atoms[neighbors[0]]));
        }
    }

    // van der waals
    std::pair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> nonbondedPair;
    foreach(nonbondedPair, interactions.nonbondedPairs()){
        addCalculation(new UffVanDerWaalsCalculation(nonbondedPair.first,
                                                     nonbondedPair.second));
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

bool UffForceField::isGroupSix(const chemkit::ForceFieldAtom *atom) const
{
    switch(atom->atom()->atomicNumber()){
        case 8:
        case 16:
        case 34:
        case 52:
        case 84:
            return true;
        default:
            return false;
    };
}
