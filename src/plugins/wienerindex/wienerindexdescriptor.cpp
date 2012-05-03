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

#include "wienerindexdescriptor.h"

#include <set>
#include <vector>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

namespace {

// Returns the graph distance between the two atoms.
int distanceBetween(const chemkit::Atom *a, const chemkit::Atom *b)
{
    int distance = 0;

    std::set<const chemkit::Atom*> visited;

    std::vector<const chemkit::Atom*> row;
    row.push_back(a);

    while(!row.empty()){
        std::vector<const chemkit::Atom*> nextRow;

        foreach(const chemkit::Atom *atom, row){
            if(atom == b){
                return distance;
            }

            visited.insert(atom);

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(visited.find(neighbor) != visited.end()){
                    continue;
                }
                else if(neighbor->isTerminalHydrogen()){
                    continue;
                }

                nextRow.push_back(neighbor);
            }
        }

        row = nextRow;
        distance++;
    }

    return 0;
}

} // end anonymous namespace

WienerIndexDescriptor::WienerIndexDescriptor()
    : chemkit::MolecularDescriptor("wiener-index")
{
    setDimensionality(2);
}

WienerIndexDescriptor::~WienerIndexDescriptor()
{
}

// Returns the wiener index for the molecule.
chemkit::Variant WienerIndexDescriptor::value(const chemkit::Molecule *molecule) const
{
    int index = 0;

    for(size_t i = 0; i < molecule->atomCount(); i++){
        const chemkit::Atom *a = molecule->atom(i);
        if(a->isTerminalHydrogen()){
            continue;
        }

        for(size_t j = i + 1; j < molecule->atomCount(); j++){
            const chemkit::Atom *b = molecule->atom(j);
            if(b->isTerminalHydrogen()){
                continue;
            }

            index += distanceBetween(a, b);
        }
    }

    return index;
}
