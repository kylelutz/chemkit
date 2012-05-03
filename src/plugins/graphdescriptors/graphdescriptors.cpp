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

#include "graphdescriptors.h"

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

                nextRow.push_back(neighbor);
            }
        }

        row = nextRow;
        distance++;
    }

    return 0;
}

} // end anonymous namespace

// === GraphDensityDescriptor ============================================== //
GraphDensityDescriptor::GraphDensityDescriptor()
    : chemkit::MolecularDescriptor("graph-density")
{
    setDimensionality(2);
}

chemkit::Variant GraphDensityDescriptor::value(const chemkit::Molecule *molecule) const
{
    if(molecule->isEmpty()){
        return chemkit::Variant();
    }

    int v = molecule->atomCount();
    int e = molecule->bondCount();

    return (2.0 * e) / (v * (v - 1.0));
}

// === GraphDiameterDescriptor ============================================= //
GraphDiameterDescriptor::GraphDiameterDescriptor()
    : chemkit::MolecularDescriptor("graph-diameter")
{
    setDimensionality(2);
}

chemkit::Variant GraphDiameterDescriptor::value(const chemkit::Molecule *molecule) const
{
    int diameter = 0;

    for(size_t i = 0; i < molecule->size(); i++){
        const chemkit::Atom *a = molecule->atom(i);

        for(size_t j = i + 1; j < molecule->size(); j++){
            const chemkit::Atom *b = molecule->atom(j);

            int distance = distanceBetween(a, b);

            if(distance > diameter){
                diameter = distance;
            }
        }
    }

    return diameter;
}

// === GraphOrderDescriptor ================================================ //
GraphOrderDescriptor::GraphOrderDescriptor()
    : chemkit::MolecularDescriptor("graph-order")
{
    setDimensionality(2);
}

chemkit::Variant GraphOrderDescriptor::value(const chemkit::Molecule *molecule) const
{
    return molecule->atomCount();
}

// === GraphRadiusDescriptor =============================================== //
GraphRadiusDescriptor::GraphRadiusDescriptor()
    : chemkit::MolecularDescriptor("graph-radius")
{
    setDimensionality(2);
}

chemkit::Variant GraphRadiusDescriptor::value(const chemkit::Molecule *molecule) const
{
    int radius = std::numeric_limits<int>::max();

    for(size_t i = 0; i < molecule->size(); i++){
        const chemkit::Atom *a = molecule->atom(i);

        int eccentricity = 0;

        for(size_t j = 0; j < molecule->size(); j++){
            const chemkit::Atom *b = molecule->atom(j);

            int distance = distanceBetween(a, b);

            if(distance > eccentricity){
                eccentricity = distance;
            }
        }

        if(eccentricity < radius){
            radius = eccentricity;
        }
    }

    return radius;
}

// === GraphSizeDescriptor ================================================= //
GraphSizeDescriptor::GraphSizeDescriptor()
    : chemkit::MolecularDescriptor("graph-size")
{
    setDimensionality(2);
}

chemkit::Variant GraphSizeDescriptor::value(const chemkit::Molecule *molecule) const
{
    return molecule->bondCount();
}
