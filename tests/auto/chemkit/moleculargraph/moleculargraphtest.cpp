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

#include "moleculargraphtest.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculargraph.h>

void MolecularGraphTest::initTestCase()
{
    m_empty = new chemkit::Molecule;
    QCOMPARE(m_empty->formula(), std::string());
    m_benzene = new chemkit::Molecule("1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
    QCOMPARE(m_benzene->formula(), std::string("C6H6"));
    m_diphenylEther = new chemkit::Molecule("1/C12H10O/c1-3-7-11(8-4-1)13-12-9-5-2-6-10-12/h1-10H", "inchi");
    QCOMPARE(m_diphenylEther->formula(), std::string("C12H10O"));
    m_ethanol = new chemkit::Molecule("1/C2H6O/c1-2-3/h3H,2H2,1H3", "inchi");
    QCOMPARE(m_ethanol->formula(), std::string("C2H6O"));
    m_octane = new chemkit::Molecule("1/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3", "inchi");
    QCOMPARE(m_octane->formula(), std::string("C8H18"));
    m_tyrosine = new chemkit::Molecule("1/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4-2-6/h1-4,8,11H,5,10H2,(H,12,13)", "inchi");
    QCOMPARE(m_tyrosine->formula(), std::string("C9H11NO3"));
}

void MolecularGraphTest::cleanupTestCase()
{
    delete m_empty;
    delete m_benzene;
    delete m_diphenylEther;
    delete m_ethanol;
    delete m_octane;
    delete m_tyrosine;
}

void MolecularGraphTest::basic()
{
    chemkit::MolecularGraph graph(m_ethanol);
    QVERIFY(graph.molecule() == m_ethanol);
    QCOMPARE(graph.atomCount(), 9U);
    QCOMPARE(graph.bondCount(), 8U);
    QCOMPARE(graph.size(), 9U);
    QCOMPARE(graph.isEmpty(), false);

    for(unsigned int i = 0; i < graph.size(); i++){
        const chemkit::Atom *atom = graph.atom(i);

        if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(graph.neighborCount(i), 4U);
        }
        else if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(graph.neighborCount(i), 1U);
        }
        else if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(graph.neighborCount(i), 2U);
        }
    }

    graph = chemkit::MolecularGraph(m_empty);
    QCOMPARE(graph.atomCount(), 0U);
    QCOMPARE(graph.bondCount(), 0U);
    QCOMPARE(graph.size(), 0U);
}

void MolecularGraphTest::cyclicGraph()
{
    chemkit::MolecularGraph *graph;

    // empty
    graph = chemkit::MolecularGraph::cyclicGraph(m_empty);
    QCOMPARE(graph->size(), 0U);
    delete graph;

    // benzene
    graph = chemkit::MolecularGraph::cyclicGraph(m_benzene);
    QCOMPARE(graph->size(), 6U);
    delete graph;

    // diphenyl ether
    graph = chemkit::MolecularGraph::cyclicGraph(m_diphenylEther);
    QCOMPARE(graph->size(), 13U);
    delete graph;

    // ethanol
    graph = chemkit::MolecularGraph::cyclicGraph(m_ethanol);
    QCOMPARE(graph->size(), 0U);
    delete graph;

    // octane
    graph = chemkit::MolecularGraph::cyclicGraph(m_octane);
    QCOMPARE(graph->size(), 0U);
    delete graph;

    // tyrosine
    graph = chemkit::MolecularGraph::cyclicGraph(m_tyrosine);
    QCOMPARE(graph->size(), 6U);
    delete graph;
}

void MolecularGraphTest::hydrogenDepletedGraph()
{
    chemkit::MolecularGraph *graph;

    // empty
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_empty);
    QCOMPARE(graph->size(), 0U);
    delete graph;

    // benzene
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_benzene);
    QCOMPARE(graph->size(), 6U);
    for(unsigned int i = 0; i < graph->size(); i++)
        QCOMPARE(graph->atom(i)->isTerminalHydrogen(), false);
    delete graph;

    // diphenyl ether
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_diphenylEther);
    QCOMPARE(graph->size(), 13U);
    for(unsigned int i = 0; i < graph->size(); i++)
        QCOMPARE(graph->atom(i)->isTerminalHydrogen(), false);
    delete graph;

    // ethanol
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_ethanol);
    QCOMPARE(graph->size(), 3U);
    for(unsigned int i = 0; i < graph->size(); i++)
        QCOMPARE(graph->atom(i)->isTerminalHydrogen(), false);
    delete graph;

    // octane
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_octane);
    QCOMPARE(graph->size(), 8U);
    for(unsigned int i = 0; i < graph->size(); i++)
        QCOMPARE(graph->atom(i)->isTerminalHydrogen(), false);
    delete graph;

    // tyrosine
    graph = chemkit::MolecularGraph::hydrogenDepletedGraph(m_tyrosine);
    QCOMPARE(graph->size(), 13U);
    for(unsigned int i = 0; i < graph->size(); i++)
        QCOMPARE(graph->atom(i)->isTerminalHydrogen(), false);
    delete graph;
}

void MolecularGraphTest::isomorphism()
{
    chemkit::MolecularGraph *benzene = chemkit::MolecularGraph::hydrogenDepletedGraph(m_benzene);
    chemkit::MolecularGraph *tyrosine = chemkit::MolecularGraph::hydrogenDepletedGraph(m_tyrosine);

    std::map<chemkit::Atom *, chemkit::Atom *> mapping = chemkit::MolecularGraph::isomorphism(benzene, tyrosine);
    QVERIFY(mapping.size() == 6);

    delete benzene;
    delete tyrosine;
}

QTEST_APPLESS_MAIN(MolecularGraphTest)
