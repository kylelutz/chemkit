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

#include "uffatomtyper.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>

// --- Construction and Destruction ---------------------------------------- //
UffAtomTyper::UffAtomTyper(const chemkit::Molecule *molecule)
    : chemkit::AtomTyper("uff")
{
    setMolecule(molecule);
}

UffAtomTyper::~UffAtomTyper()
{
}

// --- Types --------------------------------------------------------------- //
std::string UffAtomTyper::typeString(int index) const
{
    return m_types.value(index, std::string());
}

std::string UffAtomTyper::typeString(const chemkit::Atom *atom) const
{
    return typeString(atom->index());
}

// --- Internal Methods ---------------------------------------------------- //
std::string UffAtomTyper::atomType(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Hydrogen)){
        if(atom->isTerminal()){
            return "H_";
        }
        else if(atom->neighborCount() == 2){
            return "H_b";
        }
    }
    else if(atom->is(chemkit::Atom::Helium)){
        return "He4+4";
    }
    else if(atom->is(chemkit::Atom::Lithium)){
        return "Li";
    }
    else if(atom->is(chemkit::Atom::Beryllium)){
        return "Be3+2";
    }
    else if(atom->is(chemkit::Atom::Boron)){
        if(atom->neighborCount() == 2){
            return "B_2";
        }
        else if(atom->neighborCount() == 3){
            return "B_3";
        }
    }
    else if(atom->is(chemkit::Atom::Carbon)){
        if(atom->neighborCount() == 4){
            return "C_3";
        }
        else if(atom->isAromatic()){
            return "C_R";
        }
        else if(atom->neighborCount() == 3){
            return "C_2";
        }
        else if(atom->neighborCount() == 2){
            return "C_1";
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(atom->neighborCount() == 4){
            return "N_3";
        }
        else if(atom->isAromatic()){
            return "N_R";
        }
        else if(atom->neighborCount() == 3){
            return "N_2";
        }
        else if(atom->neighborCount() == 2){
            return "N_1";
        }
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(atom->isAromatic()){
            return "O_R";
        }
        else if(atom->neighborCount() == 2){
            return "O_3";
        }
        else if(atom->neighborCount() == 1){
            return "O_2";
        }
    }
    else if(atom->is(chemkit::Atom::Fluorine)){
        return "F_";
    }
    else if(atom->is(chemkit::Atom::Neon)){
        return "Ne4+4";
    }
    else if(atom->is(chemkit::Atom::Sodium)){
        return "Na";
    }
    else if(atom->is(chemkit::Atom::Magnesium)){
        return "Mg3+2";
    }
    else if(atom->is(chemkit::Atom::Aluminum)){
        return "Al3";
    }
    else if(atom->is(chemkit::Atom::Silicon)){
        return "Si3";
    }
    else if(atom->is(chemkit::Atom::Phosphorus)){
        if(atom->neighborCount() == 4){
            return "P_3+3";
        }
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(atom->neighborCount() == 4){
            return "S_3+2";
        }
        else if(atom->isAromatic()){
            return "S_R";
        }
        else if(atom->neighborCount() == 3){
            return "S_2";
        }
    }
    else if(atom->is(chemkit::Atom::Chlorine)){
        return "Cl";
    }
    else if(atom->is(chemkit::Atom::Argon)){
        return "Ar4+4";
    }
    else if(atom->is(chemkit::Atom::Potassium)){
        return "K_";
    }
    else if(atom->is(chemkit::Atom::Calcium)){
        return "Ca6+2";
    }
    else if(atom->is(chemkit::Atom::Scandium)){
        return "Sc3+3";
    }
    else if(atom->is(chemkit::Atom::Titanium)){
        return "Ti3+4";
    }
    else if(atom->is(chemkit::Atom::Vanadium)){
        return "V_3+5";
    }
    else if(atom->is(chemkit::Atom::Chromium)){
        return "Cr6+3";
    }
    else if(atom->is(chemkit::Atom::Manganese)){
        return "Mn6+2";
    }
    else if(atom->is(chemkit::Atom::Iron)){
        return "Fe3+2";
    }
    else if(atom->is(chemkit::Atom::Cobalt)){
        return "Co6+3";
    }
    else if(atom->is(chemkit::Atom::Nickel)){
        return "Ni4+2";
    }
    else if(atom->is(chemkit::Atom::Copper)){
        return "Cu3+1";
    }
    else if(atom->is(chemkit::Atom::Zinc)){
        return "Zn3+2";
    }
    else if(atom->is(chemkit::Atom::Gallium)){
        return "Ga3+3";
    }
    else if(atom->is(chemkit::Atom::Germanium)){
        return "Ge3";
    }
    else if(atom->is(chemkit::Atom::Arsenic)){
        return "As3+3";
    }
    else if(atom->is(chemkit::Atom::Selenium)){
        return "Se3+2";
    }
    else if(atom->is(chemkit::Atom::Bromine)){
        return "Br";
    }
    else if(atom->is(chemkit::Atom::Krypton)){
        return "Kr4+4";
    }
    else if(atom->is(chemkit::Atom::Rubidium)){
        return "Rb";
    }
    else if(atom->is(chemkit::Atom::Iodine)){
        return "I_";
    }

    return std::string();
}

void UffAtomTyper::assignTypes(const chemkit::Molecule *molecule)
{
    if(!molecule){
        m_types.resize(0);
        return;
    }

    m_types = QVector<std::string>(molecule->atomCount());

    for(int index = 0; index < molecule->size(); index++){
        const chemkit::Atom *atom = molecule->atom(index);

        m_types[index] = atomType(atom);
    }
}
