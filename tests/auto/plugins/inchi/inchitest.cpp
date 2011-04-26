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

#include "inchitest.h"

#include <algorithm>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/moleculefileformat.h>

void InchiTest::initTestCase()
{
    std::vector<std::string> lineFormats = chemkit::LineFormat::formats();
    QVERIFY(std::find(lineFormats.begin(), lineFormats.end(), "inchi") != lineFormats.end());
    QVERIFY(std::find(lineFormats.begin(), lineFormats.end(), "inchikey") != lineFormats.end());

    std::vector<std::string> fileFormats = chemkit::MoleculeFileFormat::formats();
    QVERIFY(std::find(fileFormats.begin(), fileFormats.end(), "inchi") != fileFormats.end());
}

void InchiTest::read()
{
    chemkit::LineFormat *inchi = chemkit::LineFormat::create("inchi");
    QVERIFY(inchi != 0);

    // empty
    chemkit::Molecule *empty = inchi->read("");
    QVERIFY(empty == 0);

    // methane
    chemkit::Molecule *methane = inchi->read("InChI=1S/CH4/h1H4");
    QVERIFY(methane != 0);
    QCOMPARE(methane->atomCount(), 5);
    QCOMPARE(methane->bondCount(), 4);
    QCOMPARE(methane->formula(), std::string("CH4"));
    delete methane;

    // ethanol
    chemkit::Molecule *ethanol = inchi->read("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
    QVERIFY(ethanol != 0);
    QCOMPARE(ethanol->atomCount(), 9);
    QCOMPARE(ethanol->bondCount(), 8);
    QCOMPARE(ethanol->formula(), std::string("C2H6O"));
    delete ethanol;

    // benzene
    chemkit::Molecule *benzene = inchi->read("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H");
    QVERIFY(benzene != 0);
    QCOMPARE(benzene->atomCount(), 12);
    QCOMPARE(benzene->bondCount(), 12);
    QCOMPARE(benzene->formula(), std::string("C6H6"));
    QCOMPARE(benzene->ringCount(), 1);
    chemkit::Ring *benzeneRing = benzene->rings()[0];
    QCOMPARE(benzeneRing->isAromatic(), true);
    delete benzene;

    delete inchi;
}

void InchiTest::write()
{
    chemkit::LineFormat *inchi = chemkit::LineFormat::create("inchi");
    QVERIFY(inchi != 0);

    chemkit::LineFormat *inchikey = chemkit::LineFormat::create("inchikey");
    QVERIFY(inchikey != 0);

    // empty molecule
    chemkit::Molecule empty;
    QCOMPARE(inchi->write(&empty), std::string(""));

    // methane molecule
    chemkit::Molecule methane;
    chemkit::Atom *methane_c1 = methane.addAtom("C");
    for(int i = 0; i < 4; i++){
        chemkit::Atom *h = methane.addAtom("H");
        methane.addBond(methane_c1, h);
    }
    QCOMPARE(methane.formula(), std::string("CH4"));
    QCOMPARE(inchi->write(&methane), std::string("InChI=1S/CH4/h1H4"));
    QCOMPARE(inchikey->write(&methane), std::string("VNWKTOKETHGBQD-UHFFFAOYSA-N"));

    // ethanol
    chemkit::Molecule ethanol;
    chemkit::Atom *ethanol_c1 = ethanol.addAtom("C");
    chemkit::Atom *ethanol_c2 = ethanol.addAtom("C");
    chemkit::Atom *ethanol_o1 = ethanol.addAtom("O");
    ethanol.addBond(ethanol_c1, ethanol_c2);
    ethanol.addBond(ethanol_c2, ethanol_o1);
    QCOMPARE(inchi->write(&ethanol), std::string("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"));
    QCOMPARE(inchikey->write(&ethanol), std::string("LFQSCWFLJHTTHZ-UHFFFAOYSA-N"));

    // benzene
    chemkit::Molecule benzene;
    chemkit::Atom *benzene_c1 = benzene.addAtom("C");
    chemkit::Atom *benzene_c2 = benzene.addAtom("C");
    chemkit::Atom *benzene_c3 = benzene.addAtom("C");
    chemkit::Atom *benzene_c4 = benzene.addAtom("C");
    chemkit::Atom *benzene_c5 = benzene.addAtom("C");
    chemkit::Atom *benzene_c6 = benzene.addAtom("C");
    benzene.addBond(benzene_c1, benzene_c2, 1);
    benzene.addBond(benzene_c2, benzene_c3, 2);
    benzene.addBond(benzene_c3, benzene_c4, 1);
    benzene.addBond(benzene_c4, benzene_c5, 2);
    benzene.addBond(benzene_c5, benzene_c6, 1);
    benzene.addBond(benzene_c6, benzene_c1, 2);
    QCOMPARE(inchi->write(&benzene), std::string("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"));
    QCOMPARE(inchikey->write(&benzene), std::string("UHOVQNZJYSORNB-UHFFFAOYSA-N"));

    delete inchi;
    delete inchikey;
}

void InchiTest::stereochemistry()
{
    chemkit::LineFormat *inchi = chemkit::LineFormat::create("inchi");
    QVERIFY(inchi != 0);

    // by default stereochemistry is on
    QCOMPARE(inchi->option("stereochemistry").toBool(), true);

    // set to false
    inchi->setOption("stereochemistry", false);
    QCOMPARE(inchi->option("stereochemistry").toBool(), false);

    // build chiral molecule
    chemkit::Molecule bromochlorofluoromethane;
    chemkit::Atom *C1 = bromochlorofluoromethane.addAtom("C");
    chemkit::Atom *Br2 = bromochlorofluoromethane.addAtom("Br");
    chemkit::Atom *Cl3 = bromochlorofluoromethane.addAtom("Cl");
    chemkit::Atom *F4 = bromochlorofluoromethane.addAtom("F");
    chemkit::Atom *H5 = bromochlorofluoromethane.addAtom("H");
    bromochlorofluoromethane.addBond(C1, Br2);
    bromochlorofluoromethane.addBond(C1, Cl3);
    bromochlorofluoromethane.addBond(C1, F4);
    bromochlorofluoromethane.addBond(C1, H5);
    QCOMPARE(inchi->write(&bromochlorofluoromethane), std::string("InChI=1S/CHBrClF/c2-1(3)4/h1H"));

    // set stererochemistry to true
    inchi->setOption("stereochemistry", true);
    QCOMPARE(inchi->option("stereochemistry").toBool(), true);

    C1->setChirality(chemkit::Atom::R);
    QCOMPARE(inchi->write(&bromochlorofluoromethane), std::string("InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1"));

    C1->setChirality(chemkit::Atom::S);
    QCOMPARE(inchi->write(&bromochlorofluoromethane), std::string("InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1"));

    C1->setChirality(chemkit::Atom::NoChirality);
    QCOMPARE(inchi->write(&bromochlorofluoromethane), std::string("InChI=1S/CHBrClF/c2-1(3)4/h1H"));

    // set stereochemistry to off
    inchi->setOption("stereochemistry", false);
    C1->setChirality(chemkit::Atom::R);
    QCOMPARE(inchi->write(&bromochlorofluoromethane), std::string("InChI=1S/CHBrClF/c2-1(3)4/h1H"));

    delete inchi;
}

void InchiTest::addHydrogens()
{
    chemkit::LineFormat *inchi = chemkit::LineFormat::create("inchi");
    QVERIFY(inchi != 0);

    // by default add-hydrogens is true
    QCOMPARE(inchi->option("add-hydrogens").toBool(), true);

    // set to false
    inchi->setOption("add-hydrogens", false);
    QCOMPARE(inchi->option("add-hydrogens").toBool(), false);

    // read octane molecule with add-hydrogens enabled
    chemkit::Molecule *octane;
    inchi->setOption("add-hydrogens", true);
    octane = inchi->read("InChI=1/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3");
    QVERIFY(octane != 0);
    QCOMPARE(octane->formula(), std::string("C8H18"));
    delete octane;

    // read octane molecule with add-hydrogens disabled
    inchi->setOption("add-hydrogens", false);
    octane = inchi->read("InChI=1/C8H18/c1-3-5-7-8-6-4-2/h3-8H2,1-2H3");
    QVERIFY(octane != 0);
    QCOMPARE(octane->formula(), std::string("C8"));
    delete octane;
}

QTEST_APPLESS_MAIN(InchiTest)
