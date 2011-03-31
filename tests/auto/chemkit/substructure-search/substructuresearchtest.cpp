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

#include <QtTest>

#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>

// The substructure-search test checks the Molecule::isSubstructureOf() method
// against a variety of small molecules. Each one is tested against every other
// to ensure that substructure mappings are correctly detected.
class SubstructureSearchTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void cleanupTestCase();

        void benzene();
        void butane();
        void cyclopropane();
        void ethane();
        void ethanol();
        void indole();
        void methane();
        void methanol();
        void phenol();
        void propane();

        void protein();

    private:
        chemkit::Molecule *m_benzene;
        chemkit::Molecule *m_butane;
        chemkit::Molecule *m_cyclopropane;
        chemkit::Molecule *m_ethane;
        chemkit::Molecule *m_ethanol;
        chemkit::Molecule *m_indole;
        chemkit::Molecule *m_methane;
        chemkit::Molecule *m_methanol;
        chemkit::Molecule *m_phenol;
        chemkit::Molecule *m_propane;
};

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
    QCOMPARE(m_benzene->isSubstructureOf(m_benzene), true);
    QCOMPARE(m_benzene->isSubstructureOf(m_butane), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_indole), true);
    QCOMPARE(m_benzene->isSubstructureOf(m_methane), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_benzene->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_benzene->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::butane()
{
    QCOMPARE(m_butane->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_butane->isSubstructureOf(m_butane), true);
    QCOMPARE(m_butane->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_butane->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_butane->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_butane->isSubstructureOf(m_indole), false);
    QCOMPARE(m_butane->isSubstructureOf(m_methane), false);
    QCOMPARE(m_butane->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_butane->isSubstructureOf(m_phenol), false);
    QCOMPARE(m_butane->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::cyclopropane()
{
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_butane), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_cyclopropane), true);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_indole), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_methane), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_phenol), false);
    QCOMPARE(m_cyclopropane->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::ethane()
{
    QCOMPARE(m_ethane->isSubstructureOf(m_benzene), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_butane), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_cyclopropane), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_ethane), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_ethanol), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_indole), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_methane), false);
    QCOMPARE(m_ethane->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_ethane->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_ethane->isSubstructureOf(m_propane), true);
}

void SubstructureSearchTest::ethanol()
{
    QCOMPARE(m_ethanol->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_butane), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_ethanol), true);
    QCOMPARE(m_ethanol->isSubstructureOf(m_indole), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_methane), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_ethanol->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_ethanol->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::indole()
{
    QCOMPARE(m_indole->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_indole->isSubstructureOf(m_butane), false);
    QCOMPARE(m_indole->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_indole->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_indole->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_indole->isSubstructureOf(m_indole), true);
    QCOMPARE(m_indole->isSubstructureOf(m_methane), false);
    QCOMPARE(m_indole->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_indole->isSubstructureOf(m_phenol), false);
    QCOMPARE(m_indole->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::methane()
{
    QCOMPARE(m_methane->isSubstructureOf(m_benzene), true);
    QCOMPARE(m_methane->isSubstructureOf(m_butane), true);
    QCOMPARE(m_methane->isSubstructureOf(m_cyclopropane), true);
    QCOMPARE(m_methane->isSubstructureOf(m_ethane), true);
    QCOMPARE(m_methane->isSubstructureOf(m_ethanol), true);
    QCOMPARE(m_methane->isSubstructureOf(m_indole), true);
    QCOMPARE(m_methane->isSubstructureOf(m_methane), true);
    QCOMPARE(m_methane->isSubstructureOf(m_methanol), true);
    QCOMPARE(m_methane->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_methane->isSubstructureOf(m_propane), true);
}

void SubstructureSearchTest::methanol()
{
    QCOMPARE(m_methanol->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_butane), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_ethanol), true);
    QCOMPARE(m_methanol->isSubstructureOf(m_indole), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_methane), false);
    QCOMPARE(m_methanol->isSubstructureOf(m_methanol), true);
    QCOMPARE(m_methanol->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_methanol->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::phenol()
{
    QCOMPARE(m_phenol->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_butane), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_cyclopropane), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_indole), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_methane), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_phenol->isSubstructureOf(m_phenol), true);
    QCOMPARE(m_phenol->isSubstructureOf(m_propane), false);
}

void SubstructureSearchTest::propane()
{
    QCOMPARE(m_propane->isSubstructureOf(m_benzene), false);
    QCOMPARE(m_propane->isSubstructureOf(m_butane), true);
    QCOMPARE(m_propane->isSubstructureOf(m_cyclopropane), true);
    QCOMPARE(m_propane->isSubstructureOf(m_ethane), false);
    QCOMPARE(m_propane->isSubstructureOf(m_ethanol), false);
    QCOMPARE(m_propane->isSubstructureOf(m_indole), true);
    QCOMPARE(m_propane->isSubstructureOf(m_methane), false);
    QCOMPARE(m_propane->isSubstructureOf(m_methanol), false);
    QCOMPARE(m_propane->isSubstructureOf(m_phenol), false);
    QCOMPARE(m_propane->isSubstructureOf(m_propane), true);
}

// This test utilizes a protein molecule from the 'alphabet.mol2' file. The
// protein consists of a single chain containing one of each of the twenty
// amino acids. Once loaded, various amino acid side chain molecules are
// created and then checked against the protein molecule to ensure that they
// are found as substructures.
void SubstructureSearchTest::protein()
{
    // read file
    chemkit::ChemicalFile file("../../../data/alphabet.mol2");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // load and verify the protein molecule
    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->atomCount(), 324);

    // indole in tryptophan
    chemkit::Molecule indole("InChI=1/C8H7N/c1-2-4-8-7(3-1)5-6-9-8/h1-6,9H", "inchi");
    QCOMPARE(indole.isSubstructureOf(molecule, chemkit::Molecule::CompareAromaticity), true);

    // phenol ring in tyrosine
    chemkit::Molecule phenol("InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H", "inchi");
    QCOMPARE(phenol.isSubstructureOf(molecule, chemkit::Molecule::CompareAromaticity), true);

    // proline ring
    chemkit::Molecule prolineRing("InChI=1/C4H9N/c1-2-4-5-3-1/h5H,1-4H2", "inchi");
    QCOMPARE(prolineRing.isSubstructureOf(molecule), true);

    // guanidine in arginine
    chemkit::Molecule guanidine("InChI=1/CH5N3/c2-1(3)4/h(H5,2,3,4)/f/h2H,3-4H2", "inchi");
    QCOMPARE(guanidine.isSubstructureOf(molecule), true);

    // lysine chain
    chemkit::Molecule lysineChain("InChI=1/C5H13N/c1-2-3-4-5-6/h2-6H2,1H3", "inchi");
    QCOMPARE(lysineChain.isSubstructureOf(molecule), true);

    // isoleucine residue
    chemkit::Molecule isoleucine("InChI=1/C6H13NO/c1-3-5(2)6(7)4-8/h4-6H,3,7H2,1-2H3", "inchi");
    QCOMPARE(isoleucine.isSubstructureOf(molecule), true);

    // cysteine side chain
    chemkit::Molecule cysteineChain("InChI=1/C2H6S/c1-2-3/h3H,2H2,1H3", "inchi");
    QCOMPARE(cysteineChain.isSubstructureOf(molecule), true);

    // methionine chain
    chemkit::Molecule methionineChain("InChI=1/C4H10S/c1-3-4-5-2/h3-4H2,1-2H3", "inchi");
    QCOMPARE(methionineChain.isSubstructureOf(molecule), true);
}

QTEST_APPLESS_MAIN(SubstructureSearchTest)
#include "substructuresearchtest.moc"
