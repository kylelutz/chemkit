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

#ifndef SMILES_H
#define SMILES_H

#include <chemkit/atom.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>

// Returns true if the atom is a member of the aromatic subset.
inline bool isAromaticElement(const chemkit::Atom *atom)
{
    switch(atom->atomicNumber()){
        case chemkit::Atom::Boron:
        case chemkit::Atom::Carbon:
        case chemkit::Atom::Nitrogen:
        case chemkit::Atom::Oxygen:
        case chemkit::Atom::Phosphorus:
        case chemkit::Atom::Sulfur:
        case chemkit::Atom::Arsenic:
        case chemkit::Atom::Selenium:
            return true;

        default:
            return false;
    }
}

inline bool isPlanarAtom(const chemkit::Atom *atom)
{
    if(atom->is(chemkit::Atom::Carbon) && atom->neighborCount() != 3){
        return false;
    }
    else if(atom->is(chemkit::Atom::Oxygen) && atom->neighborCount() != 2){
        return false;
    }
    else if(atom->is(chemkit::Atom::Sulfur) && atom->neighborCount() != 2){
        return false;
    }

    return true;
}

inline bool isAromaticRing(const chemkit::Ring *ring)
{
    foreach(const chemkit::Atom *atom, ring->atoms()){
        if(!isAromaticElement(atom)){
            return false;
        }
        else if(!isPlanarAtom(atom)){
            return false;
        }

        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(ring->contains(bond)){
                continue;
            }

            if(bond->order() == chemkit::Bond::Double && !bond->otherAtom(atom)->isInRing()){
                return false;
            }
        }
    }

    return true;
}

inline bool isAromaticAtom(const chemkit::Atom *atom)
{
    if(!isAromaticElement(atom)){
        return false;
    }

    foreach(const chemkit::Ring *ring, atom->rings()){
        if(isAromaticRing(ring)){
            return true;
        }
    }

    return false;
}

// Returns true if the element is a member of the organic subset.
inline bool isOrganicElement(const chemkit::Atom *atom)
{
    switch(atom->atomicNumber()){
        case chemkit::Atom::Boron:
        case chemkit::Atom::Carbon:
        case chemkit::Atom::Nitrogen:
        case chemkit::Atom::Oxygen:
        case chemkit::Atom::Phosphorus:
        case chemkit::Atom::Sulfur:
        case chemkit::Atom::Chlorine:
        case chemkit::Atom::Bromine:
        case chemkit::Atom::Iodine:
            return true;

        default:
            return false;
    }
}

inline bool isOrganicAtom(const chemkit::Atom *atom)
{
    return isOrganicElement(atom) && atom->formalCharge() == 0;
}

inline bool isIsotope(const chemkit::Atom *atom)
{
    if(atom->is(chemkit::Atom::Hydrogen))
        return atom->massNumber() != 1;
    else
        return atom->massNumber() != atom->atomicNumber() * 2;
}

// Returns true if the atom is an implicit hydrogen atom.
inline bool isImplicitHydrogen(const chemkit::Atom *atom)
{
    return atom->isTerminalHydrogen() &&
           atom->massNumber() == 1 &&
           !atom->isBondedTo(chemkit::Atom::Hydrogen);
}

#endif // SMILES_H
