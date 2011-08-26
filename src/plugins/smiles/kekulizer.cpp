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

#include "kekulizer.h"

#include <map>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>

#define LEMON_ONLY_TEMPLATES
#include <lemon/matching.h>
#include <lemon/list_graph.h>

namespace {

int costOfDoubleBond(const chemkit::Atom *atom)
{
    int cost = 1;

    if(atom->is(chemkit::Atom::Nitrogen)){
        cost = 2000;
    }
    else if(atom->is(chemkit::Atom::Oxygen)){
        cost = 5000;
    }
    else if(atom->is(chemkit::Atom::Sulfur)){
        cost = 4000;
    }
    else if(atom->is(chemkit::Atom::Boron)){
        cost = 100;
    }
    else if(atom->is(chemkit::Atom::Arsenic)){
        cost = 150;
    }
    else if(atom->is(chemkit::Atom::Selenium)){
        cost = 200;
    }

    cost += atom->neighborCount();

    return cost;
}

int costOfDoubleBond(const chemkit::Bond *bond)
{
    return costOfDoubleBond(bond->atom1()) + costOfDoubleBond(bond->atom2());
}

} // end anonymous namespace

// === Kekulizer =========================================================== //
void Kekulizer::kekulize(const std::vector<chemkit::Bond *> &bonds)
{
    lemon::ListGraph graph;
    lemon::ListGraph::EdgeMap<int> costs(graph);

    std::map<chemkit::Atom*, int> atomToNode;
    std::map<int, chemkit::Bond*> edgeToBond;
    foreach(chemkit::Bond *bond, bonds){
        int node1;

        std::map<chemkit::Atom*, int>::iterator iter = atomToNode.find(bond->atom1());
        if(iter != atomToNode.end()){
            node1 = iter->second;
        }
        else{
            lemon::ListGraph::Node node = graph.addNode();
            node1 = graph.id(node);
            atomToNode[bond->atom1()] = node1;
        }

        int node2;
        iter = atomToNode.find(bond->atom2());
        if(iter != atomToNode.end()){
            node2 = iter->second;
        }
        else{
            lemon::ListGraph::Node node = graph.addNode();
            node2 = graph.id(node);
            atomToNode[bond->atom2()] = node2;
        }

        lemon::ListGraph::Edge edge = graph.addEdge(graph.nodeFromId(node1), graph.nodeFromId(node2));
        int cost = costOfDoubleBond(bond);

        costs.set(edge, 10000000 - cost);

        edgeToBond[graph.id(edge)] = bond;
    }

    lemon::MaxWeightedMatching<lemon::ListGraph> matcher(graph, costs);
    matcher.run();

    for(lemon::ListGraph::EdgeIt edge(graph); edge != lemon::INVALID; ++edge){
        chemkit::Bond *bond = edgeToBond[graph.id(edge)];

        if(matcher.matching(edge)){
            bond->setOrder(chemkit::Bond::Double);
        }
        else{
            bond->setOrder(chemkit::Bond::Single);
        }
    }
}
