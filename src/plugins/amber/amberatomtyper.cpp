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

#include "amberatomtyper.h"

#include <chemkit/atom.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

// --- Construction and Destruction ---------------------------------------- //
AmberAtomTyper::AmberAtomTyper()
    : chemkit::AtomTyper("amber")
{
}

AmberAtomTyper::~AmberAtomTyper()
{
}

// --- Properties ---------------------------------------------------------- //
void AmberAtomTyper::setMolecule(const chemkit::Molecule *molecule)
{
    chemkit::AtomTyper::setMolecule(molecule);

    m_types.resize(molecule->atomCount());

    for(size_t i = 0; i < molecule->size(); i++){
        m_types[i] = assignType(molecule->atom(i));
    }
}

// --- Types --------------------------------------------------------------- //
std::string AmberAtomTyper::type(const chemkit::Atom *atom) const
{
    return m_types[atom->index()];
}

std::string AmberAtomTyper::assignType(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Hydrogen)){
        if(atom->isTerminal()){
            const chemkit::Atom *neighbor = atom->neighbor(0);

            if(neighbor->is(chemkit::Atom::Oxygen)){
                if(neighbor->neighborCount(chemkit::Atom::Hydrogen) == 2){
                    return "HW"; // hydrogen in water
                }
                else{
                    return "HO"; // hydrogen in hydroxyl
                }
            }
            else if(neighbor->is(chemkit::Atom::Carbon)){
                int electronWithdrawingGroups = 0;
                bool positiveGroup = false;
                bool aromatic = neighbor->isAromatic();
                foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                    if(secondNeighbor->formalCharge() > 0){
                        positiveGroup = true;
                    }
                    else if(secondNeighbor->is(chemkit::Atom::Oxygen) || secondNeighbor->is(chemkit::Atom::Nitrogen)){
                        electronWithdrawingGroups++;
                    }
                }

                if(positiveGroup){
                    return "HP"; // hydrogen bonded to carbon next to a positive group
                }
                else if(neighbor->neighborCount() == 2){
                    return "HZ"; // hydrogen bonded to sp carbon
                }
                else if(electronWithdrawingGroups == 0){
                    if(aromatic){
                        return "HA";
                    }
                    else{
                        return "HC";
                    }
                }
                else if(electronWithdrawingGroups == 1){
                    if(aromatic){
                        return "H4";
                    }
                    else{
                        return "H1";
                    }
                }
                else if(electronWithdrawingGroups == 2){
                    if(aromatic){
                        return "H5";
                    }
                    else{
                        return "H2";
                    }
                }
                else if(electronWithdrawingGroups == 3){
                    return "H3";
                }
            }
            else if(neighbor->is(chemkit::Atom::Nitrogen)){
                return "H"; // hydrogen bonded to nitrogen
            }
            else if(neighbor->is(chemkit::Atom::Sulfur)){
                return "HS"; // hydrogen bonded to sulfur
            }
        }
    }
    else if(atom->is(chemkit::Atom::Lithium)){
        return "Li";
    }
    else if(atom->is(chemkit::Atom::Carbon)){
        if(atom->neighborCount() == 4){
            return "CT"; // sp3 aliphatic carbon
        }
        else if(atom->neighborCount() == 3){
            if(atom->isAromatic()){
                if(atom->isInRing(5) && atom->isInRing(6)){
                    return "CB"; // sp2 aromatic carbon in 5/6 membered ring junction
                }
                else if(atom->isInRing(6)){
                    int nitrogenCount = 0;
                    foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                        if(neighbor->is(chemkit::Atom::Nitrogen) && neighbor->isInRing(6)){
                            nitrogenCount++;
                        }
                    }

                    if(nitrogenCount == 2){
                        return "CQ"; // sp2 aromatic carbon in five membered ring between two nitrogens
                    }
                    else{
                        return "CA";
                    }
                }
                else if(atom->isInRing(5) && atom->smallestRing()->atomCount(chemkit::Atom::Nitrogen) == 2){
                    return "CK"; // sp2 aromatic carbon in five membered purine ring
                }
                else{
                    return "CA"; // sp2 aromatic carbon
                }
            }
            else{
                return "C"; // sp2 carbonyl carbon
            }
        }
        else if(atom->neighborCount() == 2){
            if(atom->isBondedTo(chemkit::Atom::Nitrogen, chemkit::Bond::Triple)){
                return "CY"; // nitrile carbon
            }
            else{
                return "CZ"; // sp carbon
            }
        }
    }
    else if(atom->is(chemkit::Atom::Nitrogen)){
        if(atom->neighborCount() == 1){
            const chemkit::Atom *neighbor = atom->neighbor(0);
            if(neighbor->is(chemkit::Atom::Carbon) && atom->bondTo(neighbor)->order() == chemkit::Bond::Triple){
                return "NY"; // nitrile nitrogen
            }
        }
        else if(atom->neighborCount() == 2){
            if(atom->isAromatic() && atom->smallestRing()->size() == 5){
                return "NB"; // sp2 nitrogen in aromatic 5 membered ring with lone pair
            }
            else if(atom->isAromatic() && atom->smallestRing()->size() == 6){
                return "NC"; // sp2 nitrogen in aromatic 6 membered ring with lone pair
            }
            else{
                return "N"; // nitrogen in amide group
            }
        }
        else if(atom->neighborCount() == 3){
            if(atom->isAromatic() && atom->smallestRing()->size() == 5 && atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
                return "NA"; // sp2 nitrogen in aromatic 5 membered ring with hydrogen
            }
            else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
                return "N2"; // sp2 nitrogen in amino groups
            }
            else{
                return "N*"; // sp2 nitrogen
            }
        }
        else if(atom->neighborCount() == 4){
            if(atom->formalCharge() == 1){
                return "N3"; // nitrogen in charged amino group
            }
        }
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(atom->neighborCount() == 1){
            const chemkit::Atom *neighbor = atom->neighbor(0);
            if(neighbor->is(chemkit::Atom::Carbon)){
                bool negativeOxygen = false;
                foreach(const chemkit::Atom *secondNeighbor, neighbor->neighbors()){
                    if(secondNeighbor->is(chemkit::Atom::Oxygen) && secondNeighbor->formalCharge() < 0){
                        negativeOxygen = true;
                    }
                }

                if(negativeOxygen){
                    return "O2"; // oxygen in carboxyl and phosphate groups
                }
                else{
                    return "O"; // carbonyl oxygen
                }
            }
        }
        else if(atom->neighborCount() == 2){
            if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
                return "OW"; // sp3 oxygen in water
            }
            else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
                return "OH"; // sp3 oxygen in hydroxyl
            }
            else{
                return "OS"; // sp3 ether and ester oxygen
            }
        }
    }
    else if(atom->is(chemkit::Atom::Fluorine)){
        return "F";
    }
    else if(atom->is(chemkit::Atom::Sodium)){
        return "Na";
    }
    else if(atom->is(chemkit::Atom::Magnesium)){
        return "MG";
    }
    else if(atom->is(chemkit::Atom::Phosphorus)){
        return "P"; // phosphorus in phosphate
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
            return "SH";
        }
        else{
            return "S";
        }
    }
    else if(atom->is(chemkit::Atom::Chlorine)){
        return "Cl";
    }
    else if(atom->is(chemkit::Atom::Potassium)){
        return "K";
    }
    else if(atom->is(chemkit::Atom::Calcium)){
        return "C0";
    }
    else if(atom->is(chemkit::Atom::Iron)){
        return "FE";
    }
    else if(atom->is(chemkit::Atom::Copper)){
        return "CU";
    }
    else if(atom->is(chemkit::Atom::Bromine)){
        return "Br";
    }
    else if(atom->is(chemkit::Atom::Rubidium)){
        return "Rb";
    }
    else if(atom->is(chemkit::Atom::Iodine)){
        return "I";
    }
    else if(atom->is(chemkit::Atom::Cesium)){
        return "Cs";
    }
    else if(atom->is(chemkit::Atom::Zinc)){
        return "Zn";
    }

    return std::string();
}
