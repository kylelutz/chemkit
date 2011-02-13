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

#include "kekulizer.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>

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
void Kekulizer::kekulize(const QList<chemkit::Bond *> &bonds)
{
    lemon::ListGraph graph;
    lemon::ListGraph::EdgeMap<int> costs(graph);

    QHash<chemkit::Atom*, int> atomToNode;
    QHash<int, chemkit::Bond*> edgeToBond;
    foreach(chemkit::Bond *bond, bonds){
        int node1;
        if(atomToNode.contains(bond->atom1())){
            node1 = atomToNode[bond->atom1()];
        }
        else{
            lemon::ListGraph::Node node = graph.addNode();
            node1 = graph.id(node);
            atomToNode[bond->atom1()] = node1;
        }

        int node2;
        if(atomToNode.contains(bond->atom2())){
            node2 = atomToNode[bond->atom2()];
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
