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

#include "moleculegraphtraitstest.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculegraphtraits.h>

#include <boost/graph/bipartite.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

void MoleculeGraphTraitsTest::conceptCheck()
{
    BOOST_CONCEPT_ASSERT((boost::GraphConcept<chemkit::Molecule>));
    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<chemkit::Molecule>));
    BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<chemkit::Molecule>));
    BOOST_CONCEPT_ASSERT((boost::AdjacencyGraphConcept<chemkit::Molecule>));
    BOOST_CONCEPT_ASSERT((boost::VertexAndEdgeListGraphConcept<chemkit::Molecule>));
}

void MoleculeGraphTraitsTest::numVertices()
{
    chemkit::Molecule molecule;
    QCOMPARE(boost::num_vertices(molecule), size_t(0));

    molecule.addAtom("H");
    QCOMPARE(boost::num_vertices(molecule), size_t(1));

    molecule.addAtom("He");
    QCOMPARE(boost::num_vertices(molecule), size_t(2));

    molecule.clear();
    QCOMPARE(boost::num_vertices(molecule), size_t(0));
}

void MoleculeGraphTraitsTest::connectedComponents()
{
    chemkit::Molecule molecule("CCO", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H6O"));

    typedef std::map<const chemkit::Atom*, size_t> ComponentsMap;
    typedef boost::associative_property_map<ComponentsMap> ComponentsPropertyMap;
    ComponentsMap componentsMap;
    ComponentsPropertyMap componentsPropertyMap(componentsMap);

    size_t count = boost::connected_components(molecule, componentsPropertyMap);
    QCOMPARE(count, size_t(1));

    molecule.addAtom("He");
    count = boost::connected_components(molecule, componentsPropertyMap);
    QCOMPARE(count, size_t(2));

    molecule.addAtom("Ne");
    count = boost::connected_components(molecule, componentsPropertyMap);
    QCOMPARE(count, size_t(3));

    molecule.clear();
    count = boost::connected_components(molecule, componentsPropertyMap);
    QCOMPARE(count, size_t(0));
}

void MoleculeGraphTraitsTest::isomorphism()
{
    chemkit::Molecule ethane;
    chemkit::Atom *ethane_C1 = ethane.addAtom("C");
    chemkit::Atom *ethane_C2 = ethane.addAtom("C");
    ethane.addBond(ethane_C1, ethane_C2);

    chemkit::Molecule ethanol1;
    chemkit::Atom *ethanol1_O1 = ethanol1.addAtom("O");
    chemkit::Atom *ethanol1_C2 = ethanol1.addAtom("C");
    chemkit::Atom *ethanol1_C3 = ethanol1.addAtom("C");
    ethanol1.addBond(ethanol1_O1, ethanol1_C2);
    ethanol1.addBond(ethanol1_C2, ethanol1_C3);

    chemkit::Molecule ethanol2;
    chemkit::Atom *ethanol2_C1 = ethanol2.addAtom("C");
    chemkit::Atom *ethanol2_C2 = ethanol2.addAtom("C");
    chemkit::Atom *ethanol2_O3 = ethanol2.addAtom("O");
    ethanol2.addBond(ethanol2_C1, ethanol2_C2);
    ethanol2.addBond(ethanol2_C2, ethanol2_O3);

    QCOMPARE(boost::isomorphism(ethane, ethane), true);
    QCOMPARE(boost::isomorphism(ethane, ethanol1), false);
    QCOMPARE(boost::isomorphism(ethane, ethanol2), false);
    QCOMPARE(boost::isomorphism(ethanol1, ethanol1), true);
    QCOMPARE(boost::isomorphism(ethanol1, ethanol2), true);
    QCOMPARE(boost::isomorphism(ethanol2, ethanol2), true);
    QCOMPARE(boost::isomorphism(ethanol2, ethanol1), true);
    QCOMPARE(boost::isomorphism(ethanol1, ethane), false);
    QCOMPARE(boost::isomorphism(ethanol2, ethane), false);
}

void MoleculeGraphTraitsTest::isBipartite()
{
    chemkit::Molecule helium("[He][He]", "smiles");
    QCOMPARE(boost::is_bipartite(helium), true);

    chemkit::Molecule cyclicOzone("O1OO1", "smiles");
    QCOMPARE(boost::is_bipartite(cyclicOzone), false);

    chemkit::Molecule pentane("CCCCC", "smiles");
    QCOMPARE(boost::is_bipartite(pentane), true);

    chemkit::Molecule cyclopentane("C1CCCC1", "smiles");
    QCOMPARE(boost::is_bipartite(cyclopentane), false);

    chemkit::Molecule hexane("CCCCCC", "smiles");
    QCOMPARE(boost::is_bipartite(hexane), true);

    chemkit::Molecule cyclohexane("C1CCCCC1", "smiles");
    QCOMPARE(boost::is_bipartite(cyclohexane), true);
}

void MoleculeGraphTraitsTest::edmondsMaximumCardinalityMatching()
{
    chemkit::Molecule benzene("c1ccccc1", "smiles");
    typedef std::map<const chemkit::Atom*, const chemkit::Atom*> MatesMap;
    typedef boost::associative_property_map<MatesMap> MatesPropertyMap;
    MatesMap matesMap;
    MatesPropertyMap matesPropertyMap(matesMap);
    boost::edmonds_maximum_cardinality_matching(benzene, matesPropertyMap);
    QCOMPARE(matesMap.size(), size_t(12));
}

QTEST_APPLESS_MAIN(MoleculeGraphTraitsTest)
