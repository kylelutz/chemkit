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

#include "ambertest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>
#include <chemkit/forcefieldatom.h>

const std::string dataPath = "../../../data/";

void AmberTest::initTestCase()
{
    std::vector<std::string> forceFields = chemkit::ForceField::forceFields();
    QVERIFY(std::find(forceFields.begin(), forceFields.end(), "amber") != forceFields.end());
}

void AmberTest::adenosine()
{
    chemkit::MoleculeFile file(dataPath + "adenosine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("C10H13N5O4"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->addMolecule(molecule);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    QCOMPARE(forceField->atomCount(), 32);
    const std::vector<chemkit::ForceFieldAtom *> &atoms = forceField->atoms();
    QCOMPARE(atoms[0]->type(), std::string("CT"));
    QCOMPARE(atoms[1]->type(), std::string("OS"));
    QCOMPARE(atoms[2]->type(), std::string("CT"));
    QCOMPARE(atoms[3]->type(), std::string("CT"));
    QCOMPARE(atoms[4]->type(), std::string("OH"));
    QCOMPARE(atoms[5]->type(), std::string("CT"));
    QCOMPARE(atoms[6]->type(), std::string("OH"));
    QCOMPARE(atoms[7]->type(), std::string("CT"));
    QCOMPARE(atoms[8]->type(), std::string("OH"));
    QCOMPARE(atoms[9]->type(), std::string("N*"));
    QCOMPARE(atoms[10]->type(), std::string("CK"));
    QCOMPARE(atoms[11]->type(), std::string("NB"));
    QCOMPARE(atoms[12]->type(), std::string("CB"));
    QCOMPARE(atoms[13]->type(), std::string("CB"));
    QCOMPARE(atoms[14]->type(), std::string("NC"));
    QCOMPARE(atoms[15]->type(), std::string("CQ"));
    QCOMPARE(atoms[16]->type(), std::string("NC"));
    QCOMPARE(atoms[17]->type(), std::string("CA"));
    QCOMPARE(atoms[18]->type(), std::string("N2"));
    QCOMPARE(atoms[19]->type(), std::string("H2"));
    QCOMPARE(atoms[20]->type(), std::string("H1"));
    QCOMPARE(atoms[21]->type(), std::string("H1"));
    QCOMPARE(atoms[22]->type(), std::string("HO"));
    QCOMPARE(atoms[23]->type(), std::string("H1"));
    QCOMPARE(atoms[24]->type(), std::string("HO"));
    QCOMPARE(atoms[25]->type(), std::string("H1"));
    QCOMPARE(atoms[26]->type(), std::string("H1"));
    QCOMPARE(atoms[27]->type(), std::string("HO"));
    QCOMPARE(atoms[28]->type(), std::string("H5"));
    QCOMPARE(atoms[29]->type(), std::string("H5"));
    QCOMPARE(atoms[30]->type(), std::string("H"));
    QCOMPARE(atoms[31]->type(), std::string("H"));

    QCOMPARE(forceField->calculationCount(), 585);
    QCOMPARE(qRound(forceField->energy()), 1460);

    delete forceField;
}

void AmberTest::serine()
{
    chemkit::MoleculeFile file(dataPath + "serine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("C3H7NO3"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->addMolecule(molecule);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    QCOMPARE(forceField->atomCount(), 14);
    const std::vector<chemkit::ForceFieldAtom *> &atoms = forceField->atoms();
    QCOMPARE(atoms[0]->type(), std::string("N3"));
    QCOMPARE(atoms[1]->type(), std::string("CT"));
    QCOMPARE(atoms[2]->type(), std::string("H"));
    QCOMPARE(atoms[3]->type(), std::string("HP"));
    QCOMPARE(atoms[4]->type(), std::string("C"));
    QCOMPARE(atoms[5]->type(), std::string("CT"));
    QCOMPARE(atoms[6]->type(), std::string("H1"));
    QCOMPARE(atoms[7]->type(), std::string("H1"));
    QCOMPARE(atoms[8]->type(), std::string("OH"));
    QCOMPARE(atoms[9]->type(), std::string("HO"));
    QCOMPARE(atoms[10]->type(), std::string("O2"));
    QCOMPARE(atoms[11]->type(), std::string("O2"));
    QCOMPARE(atoms[12]->type(), std::string("H"));
    QCOMPARE(atoms[13]->type(), std::string("H"));

    QCOMPARE(forceField->calculationCount(), 118);
    QCOMPARE(qRound(forceField->energy()), 322);

    delete forceField;
}

void AmberTest::water()
{
    chemkit::Molecule water;
    chemkit::Atom *O1 = water.addAtom("O");
    chemkit::Atom *H2 = water.addAtom("H");
    chemkit::Atom *H3 = water.addAtom("H");
    water.addBond(O1, H2);
    water.addBond(O1, H3);
    O1->setPosition(1, 1, 0);
    H2->setPosition(2, 1, 0);
    H3->setPosition(1, 2, 0);

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->addMolecule(&water);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    QCOMPARE(forceField->atomCount(), 3);
    QCOMPARE(forceField->atoms()[0]->type(), std::string("OW"));
    QCOMPARE(forceField->atoms()[1]->type(), std::string("HW"));
    QCOMPARE(forceField->atoms()[2]->type(), std::string("HW"));

    QCOMPARE(forceField->calculationCount(), 3);

    QCOMPARE(qRound(forceField->energy()), 21085);

    delete forceField;
}

QTEST_APPLESS_MAIN(AmberTest)
