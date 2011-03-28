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
    foreach(const chemkit::Molecule *molecule, molecules()){
        QHash<const chemkit::Atom *, chemkit::ForceFieldAtom *> atoms;

        UffAtomTyper typer(molecule);

        foreach(const chemkit::Atom *atom, molecule->atoms()){
            chemkit::ForceFieldAtom *forceFieldAtom = new chemkit::ForceFieldAtom(this, atom);
            atoms[atom] = forceFieldAtom;
            addAtom(forceFieldAtom);
            forceFieldAtom->setType(typer.typeString(atom).c_str());
        }

        chemkit::ForceFieldInteractions interactions(molecule, this);

        // bond strech
        QPair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> bondedPair;
        foreach(bondedPair, interactions.bondedPairs()){
            addCalculation(new UffBondStrechCalculation(bondedPair.first,
                                                        bondedPair.second));
        }

        // angle bend
        QVector<const chemkit::ForceFieldAtom *> angleGroup;
        foreach(angleGroup, interactions.angleGroups()){
            addCalculation(new UffAngleBendCalculation(angleGroup[0],
                                                       angleGroup[1],
                                                       angleGroup[2]));
        }

        // torsion
        QVector<const chemkit::ForceFieldAtom *> torsionGroup;
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
                QList<chemkit::Atom *> neighbors = atom->neighbors();

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
        QPair<const chemkit::ForceFieldAtom *, const chemkit::ForceFieldAtom *> nonbondedPair;
        foreach(nonbondedPair, interactions.nonbondedPairs()){
            addCalculation(new UffVanDerWaalsCalculation(nonbondedPair.first,
                                                         nonbondedPair.second));
        }
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
