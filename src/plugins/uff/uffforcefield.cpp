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

        foreach(const chemkit::Atom *atom, molecule->atoms()){
            chemkit::ForceFieldAtom *forceFieldAtom = new chemkit::ForceFieldAtom(this, atom);
            atoms[atom] = forceFieldAtom;
            addAtom(forceFieldAtom);
            setAtomType(forceFieldAtom);
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
                QList<const chemkit::Atom *> neighbors = atom->neighbors();

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

bool UffForceField::atomsAreWithinTwoBonds(const chemkit::Atom *a, const chemkit::Atom *b) const
{
    foreach(const chemkit::Atom *neighbor, a->neighbors()){
        if(neighbor == b){
            return true;
        }
        else if(neighbor->neighbors().contains(b)){
            return true;
        }
    }

    return false;
}

void UffForceField::setAtomType(chemkit::ForceFieldAtom *forceFieldAtom)
{
    const chemkit::Atom *atom = forceFieldAtom->atom();

    switch(atom->atomicNumber()){
        case chemkit::Atom::Hydrogen:
            if(atom->isTerminal())
                forceFieldAtom->setType("H_");
            else if(atom->neighborCount() == 2)
                forceFieldAtom->setType("H_b");
            break;
        case chemkit::Atom::Helium:
            forceFieldAtom->setType("He4+4");
            break;
        case chemkit::Atom::Lithium:
            forceFieldAtom->setType("Li");
            break;
        case chemkit::Atom::Beryllium:
            forceFieldAtom->setType("Be3+2");
            break;
        case chemkit::Atom::Boron:
            if(atom->neighborCount() == 2)
                forceFieldAtom->setType("B_2");
            else if(atom->neighborCount() == 3)
                forceFieldAtom->setType("B_3");
            break;
        case chemkit::Atom::Carbon:
            if(atom->neighborCount() == 4)
                forceFieldAtom->setType("C_3");
            else if(atom->isAromatic())
                forceFieldAtom->setType("C_R");
            else if(atom->neighborCount() == 3)
                forceFieldAtom->setType("C_2");
            else if(atom->neighborCount() == 2)
                forceFieldAtom->setType("C_1");
            break;
        case chemkit::Atom::Nitrogen:
            if(atom->neighborCount() == 4)
                forceFieldAtom->setType("N_3");
            else if(atom->isAromatic())
                forceFieldAtom->setType("N_R");
            else if(atom->neighborCount() == 3)
                forceFieldAtom->setType("N_2");
            else if(atom->neighborCount() == 2)
                forceFieldAtom->setType("N_1");
            break;
        case chemkit::Atom::Oxygen:
            if(atom->isAromatic())
                forceFieldAtom->setType("O_R");
            else if(atom->neighborCount() == 2)
                forceFieldAtom->setType("O_3");
            else if(atom->neighborCount() == 1)
                forceFieldAtom->setType("O_2");
            break;
        case chemkit::Atom::Fluorine:
            forceFieldAtom->setType("F_");
            break;
        case chemkit::Atom::Neon:
            forceFieldAtom->setType("Ne4+4");
            break;
        case chemkit::Atom::Sodium:
            forceFieldAtom->setType("Na");
            break;
        case chemkit::Atom::Magnesium:
            forceFieldAtom->setType("Mg3+2");
            break;
        case chemkit::Atom::Aluminum:
            forceFieldAtom->setType("Al3");
            break;
        case chemkit::Atom::Silicon:
            forceFieldAtom->setType("Si3");
            break;
        case chemkit::Atom::Phosphorus:
            if(atom->neighborCount() == 4)
                forceFieldAtom->setType("P_3+3");
            break;
        case chemkit::Atom::Sulfur:
            if(atom->neighborCount() == 4)
                forceFieldAtom->setType("S_3+2");
            else if(atom->isAromatic())
                forceFieldAtom->setType("S_R");
            else if(atom->neighborCount() == 3)
                forceFieldAtom->setType("S_2");
            break;
        case chemkit::Atom::Chlorine:
            forceFieldAtom->setType("Cl");
            break;
        case chemkit::Atom::Argon:
            forceFieldAtom->setType("Ar4+4");
            break;
        case chemkit::Atom::Potassium:
            forceFieldAtom->setType("K_");
            break;
        case chemkit::Atom::Calcium:
            forceFieldAtom->setType("Ca6+2");
            break;
        case chemkit::Atom::Scandium:
            forceFieldAtom->setType("Sc3+3");
            break;
        case chemkit::Atom::Titanium:
            forceFieldAtom->setType("Ti3+4");
            break;
        case chemkit::Atom::Vanadium:
            forceFieldAtom->setType("V_3+5");
            break;
        case chemkit::Atom::Chromium:
            forceFieldAtom->setType("Cr6+3");
            break;
        case chemkit::Atom::Manganese:
            forceFieldAtom->setType("Mn6+2");
            break;
        case chemkit::Atom::Iron:
            forceFieldAtom->setType("Fe3+2");
            break;
        case chemkit::Atom::Cobalt:
            forceFieldAtom->setType("Co6+3");
            break;
        case chemkit::Atom::Nickel:
            forceFieldAtom->setType("Ni4+2");
            break;
        case chemkit::Atom::Copper:
            forceFieldAtom->setType("Cu3+1");
            break;
        case chemkit::Atom::Zinc:
            forceFieldAtom->setType("Zn3+2");
            break;
        case chemkit::Atom::Gallium:
            forceFieldAtom->setType("Ga3+3");
            break;
        case chemkit::Atom::Germanium:
            forceFieldAtom->setType("Ge3");
            break;
        case chemkit::Atom::Arsenic:
            forceFieldAtom->setType("As3+3");
            break;
        case chemkit::Atom::Selenium:
            forceFieldAtom->setType("Se3+2");
            break;
        case chemkit::Atom::Bromine:
            forceFieldAtom->setType("Br");
            break;
        case chemkit::Atom::Krypton:
            forceFieldAtom->setType("Kr4+4");
            break;
        case chemkit::Atom::Rubidium:
            forceFieldAtom->setType("Rb");
            break;
        case chemkit::Atom::Iodine:
            forceFieldAtom->setType("I_");
            break;

        default:
            break;
    };
}
