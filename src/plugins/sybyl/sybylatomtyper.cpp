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

#include "sybylatomtyper.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>

// --- Construction and Destruction ---------------------------------------- //
SybylAtomTyper::SybylAtomTyper()
    : chemkit::AtomTyper("sybyl")
{
}

SybylAtomTyper::~SybylAtomTyper()
{
}

// --- Types --------------------------------------------------------------- //
std::string SybylAtomTyper::typeString(int index) const
{
    return m_types.value(index);
}

std::string SybylAtomTyper::typeString(const chemkit::Atom *atom) const
{
    return typeString(atom->index());
}

// --- Internal Methods ---------------------------------------------------- //
std::string SybylAtomTyper::atomType(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Hydrogen)){
        return "H"; // hydrogen
    }
    else if(atom->is(chemkit::Atom::Lithium)){
        return "Li"; // lithium
    }
    else if(atom->is(chemkit::Atom::Carbon)){
        if(atom->isAromatic()){
            return "C.ar"; // aromatic carbon
        }
        else if(atom->neighborCount() == 2){
            return "C.1"; // sp carbon
        }
        else if(atom->neighborCount() == 3){
            return "C.2"; // sp2 carbon
        }
        else if(atom->neighborCount() == 4){
            return "C.3"; // sp3 carbon
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(atom->isAromatic()){
            return "N.ar"; // aromatic nitrogen
        }
        else if(atom->neighborCount() == 1){
            return "N.1"; // sp nitrogen
        }
        else if(atom->neighborCount() == 2){
            return "N.2"; // sp2 nitrogen
        }
        else if(atom->neighborCount() == 3){
            return "N.3"; // sp3 nitrogen
        }
        else if(atom->neighborCount() == 4){
            return "N.4"; // sp3 positively charged nitrogen
        }
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(atom->neighborCount() == 1){
            return "O.2"; // sp2 oxygen
        }
        else if(atom->neighborCount() == 2){
            return "O.3"; // sp3 oxygen
        }
    }
    else if(atom->is(chemkit::Atom::Fluorine)){
        return "F"; // fluorine
    }
    else if(atom->is(chemkit::Atom::Sodium)){
        return "Na"; // sodium
    }
    else if(atom->is(chemkit::Atom::Magnesium)){
        return "Mg"; // magnesium
    }
    else if(atom->is(chemkit::Atom::Aluminum)){
        return "Al"; // aluminum
    }
    else if(atom->is(chemkit::Atom::Silicon)){
        return "Si"; // silicon
    }
    else if(atom->is(chemkit::Atom::Phosphorus)){
        if(atom->neighborCount() == 3){
            return "P.3"; // sp3 phosphorus
        }
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(atom->neighborCount() == 3){
            int doubleBondedOxygens = 0;

            foreach(const chemkit::Bond *bond, atom->bonds()){
                const chemkit::Atom *neighbor = bond->otherAtom(atom);

                if(neighbor->is(chemkit::Atom::Oxygen) && bond->order() == chemkit::Bond::Double){
                    doubleBondedOxygens++;
                }
            }

            if(doubleBondedOxygens == 1){
                return "S.o"; // sulfoxide sulfur
            }
            else if(doubleBondedOxygens == 2){
                return "S.o2"; // sulfone sulfur
            }
        }
        else if(atom->neighborCount() == 1){
            return "S.2"; // sp2 sulfur
        }
        else if(atom->neighborCount() == 2){
            return "S.3"; // sp3 sulfur
        }
    }
    else if(atom->is(chemkit::Atom::Chlorine)){
        return "Cl"; // chlorine
    }
    else if(atom->is(chemkit::Atom::Potassium)){
        return "K"; // potassium
    }
    else if(atom->is(chemkit::Atom::Calcium)){
        return "Ca"; // calcium
    }
    else if(atom->is(chemkit::Atom::Chromium)){
        if(atom->isBondedTo(chemkit::Atom::Oxygen)){
            return "Cr.oh"; // hydroxy chromium
        }
        else{
            return "Cr.th"; // chromium
        }
    }
    else if(atom->is(chemkit::Atom::Manganese)){
        return "Mn"; // manganese
    }
    else if(atom->is(chemkit::Atom::Iron)){
        return "Fe"; // iron
    }
    else if(atom->is(chemkit::Atom::Cobalt)){
        if(atom->isBondedTo(chemkit::Atom::Oxygen)){
            return "Co.oh"; // hydroxy cobalt
        }
    }
    else if(atom->is(chemkit::Atom::Copper)){
        return "Cu"; // copper
    }
    else if(atom->is(chemkit::Atom::Zinc)){
        return "Zn"; // zinc
    }
    else if(atom->is(chemkit::Atom::Selenium)){
        return "Se"; // selenium
    }
    else if(atom->is(chemkit::Atom::Bromine)){
        return "Br"; // bromine
    }
    else if(atom->is(chemkit::Atom::Molybdenum)){
        return "Mo"; // molybdenum
    }
    else if(atom->is(chemkit::Atom::Tin)){
        return "Sn"; // tin
    }
    else if(atom->is(chemkit::Atom::Iodine)){
        return "I"; // iodine
    }

    return std::string();
}

void SybylAtomTyper::assignTypes(const chemkit::Molecule *molecule)
{
    if(!molecule){
        m_types.resize(0);
        return;
    }

    m_types = QVector<std::string>(molecule->atomCount());

    for(int i = 0; i < molecule->atomCount(); i++){
        m_types[i] = atomType(molecule->atom(i));
    }
}
