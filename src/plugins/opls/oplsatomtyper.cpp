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

// --- Properties ---------------------------------------------------------- //
void OplsAtomTyper::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::AtomTyper::setMolecule(molecule);

    if(!molecule){
        m_typeNumbers.resize(0);
        return;
    }

    m_typeNumbers = std::vector<int>(molecule->atomCount());

    for(size_t index = 0; index < molecule->atomCount(); index++){
        const chemkit::Atom *atom = molecule->atom(index);

        // hydrogen
        if(atom->is(chemkit::Atom::Hydrogen)){
            if(atom->isTerminal()){
                const chemkit::Atom *neighbor = atom->neighbor(0);

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
                const chemkit::Atom *neighbor = atom->neighbor(0);
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

// --- Types --------------------------------------------------------------- //
std::string OplsAtomTyper::type(const chemkit::Atom *atom) const
{
    return boost::lexical_cast<std::string>(m_typeNumbers[atom->index()]);
}

void OplsAtomTyper::setTypeNumber(int index, int typeNumber)
{
    m_typeNumbers[index] = typeNumber;
}

int OplsAtomTyper::typeNumber(const chemkit::Atom *atom) const
{
    return m_typeNumbers[atom->index()];
}
