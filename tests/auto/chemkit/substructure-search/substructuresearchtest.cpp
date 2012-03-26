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

// The substructure-search test checks the SubstructureQuery::matches() method
// against a variety of small molecules. Each one is tested against every other
// to ensure that substructure mappings are correctly detected.

#include "substructuresearchtest.h"

#include <boost/make_shared.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/substructurequery.h>

void SubstructureSearchTest::initTestCase()
{
    m_benzene = new chemkit::Molecule("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
    QCOMPARE(m_benzene->formula(), std::string("C6H6"));
    m_butane = new chemkit::Molecule("InChI=1/C4H10/c1-3-4-2/h3-4H2,1-2H3", "inchi");
    QCOMPARE(m_butane->formula(), std::string("C4H10"));
    m_cyclopropane = new chemkit::Molecule("InChI=1/C3H6/c1-2-3-1/h1-3H2", "inchi");
    QCOMPARE(m_cyclopropane->formula(), std::string("C3H6"));
    m_ethane = new chemkit::Molecule("InChI=1/C2H6/c1-2/h1-2H3", "inchi");
    QCOMPARE(m_ethane->formula(), std::string("C2H6"));
    m_ethanol = new chemkit::Molecule("InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3", "inchi");
    QCOMPARE(m_ethanol->formula(), std::string("C2H6O"));
    m_indole = new chemkit::Molecule("InChI=1/C8H7N/c1-2-4-8-7(3-1)5-6-9-8/h1-6,9H", "inchi");
    QCOMPARE(m_indole->formula(), std::string("C8H7N"));
    m_methane = new chemkit::Molecule("InChI=1/CH4/h1H4", "inchi");
    QCOMPARE(m_methane->formula(), std::string("CH4"));
    m_methanol = new chemkit::Molecule("InChI=1/CH4O/c1-2/h2H,1H3", "inchi");
    QCOMPARE(m_methanol->formula(), std::string("CH4O"));
    m_propane = new chemkit::Molecule("InChI=1/C3H8/c1-3-2/h3H2,1-2H3", "inchi");
    QCOMPARE(m_propane->formula(), std::string("C3H8"));
    m_phenol = new chemkit::Molecule("InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H", "inchi");
    QCOMPARE(m_phenol->formula(), std::string("C6H6O"));
}

void SubstructureSearchTest::cleanupTestCase()
{
    delete m_benzene;
    delete m_butane;
    delete m_cyclopropane;
    delete m_ethane;
    delete m_ethanol;
    delete m_indole;
    delete m_methane;
    delete m_methanol;
    delete m_phenol;
    delete m_propane;
}

void SubstructureSearchTest::benzene()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_benzene));
    QCOMPARE(query.matches(m_benzene), true);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), true);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::butane()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_butane));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), true);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), false);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), false);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::cyclopropane()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_cyclopropane));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), true);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), false);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), false);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::ethane()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_ethane));
    QCOMPARE(query.matches(m_benzene), true);
    QCOMPARE(query.matches(m_butane), true);
    QCOMPARE(query.matches(m_cyclopropane), true);
    QCOMPARE(query.matches(m_ethane), true);
    QCOMPARE(query.matches(m_ethanol), true);
    QCOMPARE(query.matches(m_indole), true);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), true);
}

void SubstructureSearchTest::ethanol()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_ethanol));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), true);
    QCOMPARE(query.matches(m_indole), false);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::indole()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_indole));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), true);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), false);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::methane()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_methane));
    QCOMPARE(query.matches(m_benzene), true);
    QCOMPARE(query.matches(m_butane), true);
    QCOMPARE(query.matches(m_cyclopropane), true);
    QCOMPARE(query.matches(m_ethane), true);
    QCOMPARE(query.matches(m_ethanol), true);
    QCOMPARE(query.matches(m_indole), true);
    QCOMPARE(query.matches(m_methane), true);
    QCOMPARE(query.matches(m_methanol), true);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), true);
}

void SubstructureSearchTest::methanol()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_methanol));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), true);
    QCOMPARE(query.matches(m_indole), false);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), true);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::phenol()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_phenol));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), false);
    QCOMPARE(query.matches(m_cyclopropane), false);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), false);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), true);
    QCOMPARE(query.matches(m_propane), false);
}

void SubstructureSearchTest::propane()
{
    chemkit::SubstructureQuery query(boost::make_shared<chemkit::Molecule>(*m_propane));
    QCOMPARE(query.matches(m_benzene), false);
    QCOMPARE(query.matches(m_butane), true);
    QCOMPARE(query.matches(m_cyclopropane), true);
    QCOMPARE(query.matches(m_ethane), false);
    QCOMPARE(query.matches(m_ethanol), false);
    QCOMPARE(query.matches(m_indole), true);
    QCOMPARE(query.matches(m_methane), false);
    QCOMPARE(query.matches(m_methanol), false);
    QCOMPARE(query.matches(m_phenol), false);
    QCOMPARE(query.matches(m_propane), true);
}

// This test utilizes a protein molecule from the 'alphabet.mol2' file. The
// protein consists of a single chain containing one of each of the twenty
// amino acids. Once loaded, various amino acid side chain molecules are
// created and then checked against the protein molecule to ensure that they
// are found as substructures.
void SubstructureSearchTest::protein()
{
    // read file
    chemkit::MoleculeFile file("../../../data/alphabet.mol2");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // load and verify the protein molecule
    const boost::shared_ptr<chemkit::Molecule> &molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->atomCount(), size_t(324));

    // create query object
    chemkit::SubstructureQuery query;
    query.setFlags(chemkit::SubstructureQuery::CompareAromaticity);

    // indole in tryptophan
    query.setMolecule("InChI=1/C8H7N/c1-2-4-8-7(3-1)5-6-9-8/h1-6,9H", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // phenol ring in tyrosine
    query.setMolecule("InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // proline ring
    query.setMolecule("InChI=1/C4H9N/c1-2-4-5-3-1/h5H,1-4H2", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // guanidine in arginine
    query.setMolecule("InChI=1/CH5N3/c2-1(3)4/h(H5,2,3,4)/f/h2H,3-4H2", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // lysine chain
    query.setMolecule("InChI=1/C5H13N/c1-2-3-4-5-6/h2-6H2,1H3", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // isoleucine residue
    query.setMolecule("InChI=1/C6H13NO/c1-3-5(2)6(7)4-8/h4-6H,3,7H2,1-2H3", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // cysteine side chain
    query.setMolecule("InChI=1/C2H6S/c1-2-3/h3H,2H2,1H3", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);

    // methionine chain
    query.setMolecule("InChI=1/C4H10S/c1-3-4-5-2/h3-4H2,1-2H3", "inchi");
    QCOMPARE(query.matches(molecule.get()), true);
}

QTEST_APPLESS_MAIN(SubstructureSearchTest)
