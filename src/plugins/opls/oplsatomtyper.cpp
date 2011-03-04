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

#include "oplsatomtyper.h"

// --- Construction and Destruction ---------------------------------------- //
OplsAtomTyper::OplsAtomTyper(const chemkit::Molecule *molecule)
    : chemkit::AtomTyper("opls")
{
    setMolecule(molecule);
}

OplsAtomTyper::~OplsAtomTyper()
{
}

// --- Types --------------------------------------------------------------- //
void OplsAtomTyper::setTypeNumber(int index, int typeNumber)
{
    m_typeNumbers[index] = typeNumber;
}

int OplsAtomTyper::typeNumber(const chemkit::Atom *atom) const
{
    return m_typeNumbers[atom->index()];
}

QString OplsAtomTyper::typeString(const chemkit::Atom *atom) const
{
    return QString::number(typeNumber(atom));
}

void OplsAtomTyper::assignTypes(const chemkit::Molecule *molecule)
{
    if(!molecule){
        m_typeNumbers.resize(0);
        return;
    }

    m_typeNumbers = QVector<int>(molecule->atomCount());

    for(int index = 0; index < molecule->atomCount(); index++){
        const chemkit::Atom *atom = molecule->atom(index);

        // hydrogen
        if(atom->is(chemkit::Atom::Hydrogen)){
            if(atom->isTerminal()){
                const chemkit::Atom *neighbor = atom->neighbors()[0];

                if(neighbor->is(chemkit::Atom::Oxygen)){
                    if(neighbor->neighborCount() == 2 &&
                       neighbor->neighborCount(chemkit::Atom::Hydrogen) == 2){
                        setTypeNumber(index, 76); // SPC hydrogen in water (HW)
                    }
                    else{
                        setTypeNumber(index, 94); // hydrogen in alcohol (HO)
                    }
                }
                else if(neighbor->is(chemkit::Atom::Carbon)){
                    setTypeNumber(index, 82); // alkane C-H
                }
                else if(neighbor->is(chemkit::Atom::Nitrogen)){
                    if(neighbor->neighborCount(chemkit::Atom::Hydrogen) == 3){
                        setTypeNumber(index, 70); // hydrogen in ammonia (H)
                    }
                }
            }
        }
        // helium
        else if(atom->is(chemkit::Atom::Helium)){
            setTypeNumber(index, 43); // helium atom
        }
        // lithium
        else if(atom->is(chemkit::Atom::Lithium)){
            setTypeNumber(index, 345); // lithium 1+ ion (Li)
        }
        // carbon
        else if(atom->is(chemkit::Atom::Carbon)){
            if(atom->neighborCount() == 4){
                if(atom->neighborCount(chemkit::Atom::Carbon) == 2){
                    setTypeNumber(index, 78); // alkane -CH2-
                }
                else if(atom->neighborCount(chemkit::Atom::Carbon) == 1 &&
                        atom->neighborCount(chemkit::Atom::Hydrogen) == 3){
                    setTypeNumber(index, 77); // alkane -CH3
                }
                else if(atom->neighborCount(chemkit::Atom::Oxygen) == 1){
                    setTypeNumber(index, 96); // alcohol CH3OH
                }
            }
            else if(atom->neighborCount() == 3){
                if(atom->isAromatic()){
                    setTypeNumber(index, 87); // aromatic carbon
                }
            }
        }
        // nitrogen
        else if(atom->is(chemkit::Atom::Nitrogen)){
            if(atom->neighborCount() == 3){
                if(atom->neighborCount(chemkit::Atom::Hydrogen) == 3){
                    setTypeNumber(index, 69); // nitrogen in ammonia (NT)
                }
            }
        }
        // oxygen
        else if(atom->is(chemkit::Atom::Oxygen)){
            if(atom->neighborCount() == 1){
                const chemkit::Atom *neighbor = atom->neighbors()[0];
                const chemkit::Bond *neighborBond = atom->bonds()[0];

                if(neighbor->is(chemkit::Atom::Carbon) && neighborBond->order() == chemkit::Bond::Double){
                    setTypeNumber(index, 220); // ketone C=O (O)
                }
            }
            else if(atom->neighborCount() == 2){
                if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
                    setTypeNumber(index, 75); // SPC oxygen in water (OW)
                }
                else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
                    setTypeNumber(index, 93); // oxygen in alcohol (OH)
                }
            }
        }
        // fluorine
        else if(atom->is(chemkit::Atom::Fluorine)){
            if(atom->formalCharge() < 0){
                setTypeNumber(index, 340); // fluoride ion (F)
            }
        }
        // neon
        else if(atom->is(chemkit::Atom::Neon)){
            setTypeNumber(index, 44); // neon atom
        }
        // sodium
        else if(atom->is(chemkit::Atom::Sodium)){
            setTypeNumber(index, 346); // sodium ion
        }
        // magnesium
        else if(atom->is(chemkit::Atom::Magnesium)){
            setTypeNumber(index, 350); // magnesium ion (Mg)
        }
        // phosphorus
        else if(atom->is(chemkit::Atom::Phosphorus)){
            if(atom->neighborCount() == 4){
                if(atom->neighborCount(chemkit::Atom::Oxygen) > 0){
                    setTypeNumber(index, 378); // phosphate P
                }
            }
        }
        // sulfur
        else if(atom->is(chemkit::Atom::Sulfur)){
            if(atom->neighborCount() == 2){
                if(atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
                    setTypeNumber(index, 139); // sulfur in thiol (SH)
                }
                else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
                    setTypeNumber(index, 140); // sulfur in hydrogen sulfide (SH)
                }
                else if(atom->neighborCount(chemkit::Atom::Sulfur) == 1){
                    setTypeNumber(index, 142); // disulfide -S-S- (S)
                }
                else{
                    setTypeNumber(index, 141); // sulfide -S- (S)
                }
            }
        }
        // chlorine
        else if(atom->is(chemkit::Atom::Chlorine)){
            if(atom->formalCharge() < 0){
                setTypeNumber(index, 341); // chloride ion (Cl)
            }
        }
        // argon
        else if(atom->is(chemkit::Atom::Argon)){
            setTypeNumber(index, 45); // argon atom
        }
        // potassium
        else if(atom->is(chemkit::Atom::Potassium)){
            setTypeNumber(index, 347); // potassium 1+ ion (K)
        }
        // calcium
        else if(atom->is(chemkit::Atom::Calcium)){
            setTypeNumber(index, 351); // calcium 2+ ion (Ca)
        }
        // zinc
        else if(atom->is(chemkit::Atom::Zinc)){
            if(atom->formalCharge() == 2){
                setTypeNumber(index, 834); // zinc 2+ ion (Zn)
            }
        }
        // bromine
        else if(atom->is(chemkit::Atom::Bromine)){
            if(atom->formalCharge() < 0){
                setTypeNumber(index, 342); // bromide ion (Br)
            }
        }
        // krypton
        else if(atom->is(chemkit::Atom::Krypton)){
            setTypeNumber(index, 46); // krypton atom
        }
        // iodine
        else if(atom->is(chemkit::Atom::Iodine)){
            setTypeNumber(index, 343); // iodide ion (I)
        }
        // xenon
        else if(atom->is(chemkit::Atom::Xenon)){
            setTypeNumber(index, 47); // xenon atom
        }
    }
}
