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

#include "mmffaromaticitymodel.h"

#include <chemkit/atom.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>

MmffAromaticityModel::MmffAromaticityModel()
    : chemkit::AromaticityModel("mmff")
{
}

MmffAromaticityModel::~MmffAromaticityModel()
{
}

bool MmffAromaticityModel::isAromaticRing(const chemkit::Ring *ring) const
{
    if(ring->size() != 5 && ring->size() != 6)
        return false;

    int piCount = piElectronCount(ring);

    // count exocyclic aromatic bonds
    foreach(const chemkit::Atom *atom, ring->atoms()){
        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(ring->contains(bond))
                continue;

            if(bond->order() == chemkit::Bond::Double){
                foreach(const chemkit::Ring *otherRing, bond->rings()){
                    if(otherRing == ring)
                        continue;

                    else if(piElectronCount(otherRing) == 6)
                        piCount += 1;
                }
            }
        }
    }

    return piCount == 6;
}

int MmffAromaticityModel::piElectronCount(const chemkit::Ring *ring) const
{
    int piElectronCount = 0;

    // ring lone pair donors
    foreach(const chemkit::Atom *atom, ring->atoms()){
        if(ring->size() == 5){
            if(atom->is(chemkit::Atom::Nitrogen) &&
               atom->neighborCount() == 3 &&
               atom->valence() == 3){
                piElectronCount += 2;
                break;
            }
            else if(atom->is(chemkit::Atom::Nitrogen) &&
                    atom->neighborCount() == 2 &&
                    atom->valence() == 2){
                piElectronCount += 2;
                break;
            }
            else if((atom->is(chemkit::Atom::Oxygen) || atom->is(chemkit::Atom::Sulfur)) &&
                    atom->neighborCount() == 2){
                piElectronCount += 2;
                break;
            }
        }
    }

    // ring double bonds
    foreach(const chemkit::Bond *bond, ring->bonds()){
        if(bond->order() == chemkit::Bond::Double){
            piElectronCount += 2;
        }
    }

    return piElectronCount;
}
