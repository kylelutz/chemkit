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

#include "uffatomtyper.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
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

// --- Properties ---------------------------------------------------------- //
void UffAtomTyper::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::AtomTyper::setMolecule(molecule);

    if(!molecule){
        m_types.resize(0);
        return;
    }

    m_types = std::vector<std::string>(molecule->atomCount());

    for(size_t index = 0; index < molecule->size(); index++){
        const chemkit::Atom *atom = molecule->atom(index);

        m_types[index] = atomType(atom);
    }
}

// --- Types --------------------------------------------------------------- //
std::string UffAtomTyper::type(const chemkit::Atom *atom) const
{
    return m_types[atom->index()];
}

// --- Interaction Types --------------------------------------------------- //
int UffAtomTyper::bondedInteractionType(const chemkit::Atom *a,
                                        const chemkit::Atom *b) const
{
    const chemkit::Bond *bond = a->bondTo(b);

    std::string typeA = type(a);
    std::string typeB = type(b);

    if((typeA.length() > 2 && typeA[2] == 'R') && (typeB.length() > 2 && typeB[2] == 'R')){
        return Resonant;
    }
    else{
        return bond->order();
    }
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
        else if(atom->neighborCount() == 1 &&
                atom->valence() == 3){
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
