/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "tpsadescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

TpsaDescriptor::TpsaDescriptor()
    : chemkit::MolecularDescriptor("tpsa")
{
    setDimensionality(1);
}

// Returns the polar surface area contribution for the atom. The
// contribution values are listed in table 1 in [Ertl 2000].
chemkit::Real TpsaDescriptor::polarSurfaceAreaContribution(const chemkit::Atom *atom) const
{
    if(atom->is(chemkit::Atom::Nitrogen)){
        if(atom->isAromatic()){
            if(atom->isBondedTo(chemkit::Atom::Hydrogen)){
                if(atom->formalCharge() == 1){
                    return 14.14; // [nH+](:*):*
                }
                else{
                    return 15.79; // [nH](:*):*
                }
            }
            else if(atom->neighborCount() == 3){
                return 4.93; // [n](-*)(:*):*
            }
            else{
                return 12.89; // [n](:*):*
            }
        }
        if(atom->isInRing() &&
           atom->smallestRing()->size() == 3){
            if(atom->isBondedTo(chemkit::Atom::Hydrogen)){
                return 21.94; // [NH]1-*-*-1
            }
            else{
                return 3.01; // [N]1(-*)-*-*-1
            }
        }

        // check for triple bond
        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(bond->is(chemkit::Bond::Triple)){
                return 23.79; // [N]#
            }
        }

        if(atom->neighborCount(chemkit::Atom::Hydrogen) == 3 &&
           atom->formalCharge() == 1){
            return 27.64; // [NH3+]-*
        }
        else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 2){
            if(atom->formalCharge() == 1){
                // check for double bond
                foreach(const chemkit::Bond *bond, atom->bonds()){
                    if(bond->order() == chemkit::Bond::Double){
                        return 25.59; // [NH2+]=*
                    }
                }

                return 16.61; // [NH2+](-*)-*
            }
            else{
                return 26.02; // [NH2]-*
            }
        }
        else if(atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
            return 12.03; // [NH](-*)-*
        }

        // check for nitro group
        if(atom->formalCharge() == +1 &&
           atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Single) &&
           atom->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double)){
            return 11.68; // [N](-*)(=*)=*
        }

        // check for double bond
        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(bond->is(chemkit::Bond::Double)){
                return 12.36; // [N](-*)=*
            }
        }

        return 3.24; // [N](-*)(-*)-*
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        if(atom->isAromatic()){
            return 13.14; // [o](:*):*
        }
        else if(atom->isInRing() &&
                atom->smallestRing()->size() == 3){
            return 12.53; // [O]1-*-*-1
        }
        else if(atom->isTerminal() &&
                atom->bond(0)->is(chemkit::Bond::Double)){
            return 17.07; // [O]=*
        }
        else if(atom->isBondedTo(chemkit::Atom::Hydrogen)){
            return 20.23; // [OH]-*
        }
        else if(atom->isTerminal()){
            const chemkit::Atom *neighbor = atom->neighbor(0);
            if(neighbor->is(chemkit::Atom::Nitrogen) &&
               neighbor->isBondedTo(chemkit::Atom::Oxygen, chemkit::Bond::Double) &&
               atom->formalCharge() == -1){
                // nitro group
                return 17.07; // [O]=*
            }

            return 23.06; // [O-]-*
        }
        else{
            return 9.23; // [O](-*)-*
        }
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        if(atom->isAromatic()){
            // check for double bond
            foreach(const chemkit::Bond *bond, atom->bonds()){
                if(bond->is(chemkit::Bond::Double)){
                    return 21.70; // [s](=*)(:*):*
                }
            }

            return 28.24; // [s](:*):*
        }
        else if(atom->isTerminal() &&
                atom->bond(0)->is(chemkit::Bond::Double)){
            return 32.09; // [S]=*
        }
        else if(atom->neighborCount() == 2 &&
                atom->isBondedTo(chemkit::Atom::Hydrogen)){
            return 38.80; // [SH]-*
        }
    }
    else if(atom->is(chemkit::Atom::Phosphorus)){
        bool hasDoubleBond = false;
        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(bond->is(chemkit::Bond::Double)){
                hasDoubleBond = true;
                break;
            }
        }

        if(hasDoubleBond){
            if(atom->neighborCount() == 2){
                return 34.14; // [P](-*)=*
            }
            else if(atom->neighborCount() == 4 &&
                    atom->isBondedTo(chemkit::Atom::Hydrogen)){
                return 23.47; // [PH](-*)(-*)=*
            }
            else if(atom->neighborCount() == 4){
                return 9.81; // [P](-*)(-*)(-*)=*
            }
        }
        else{
            return 13.59; // [P](-*)(-*)-*
        }
    }

    return 0;
}

chemkit::Variant TpsaDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real tpsa = 0;

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        tpsa += polarSurfaceAreaContribution(atom);
    }

    return tpsa;
}
