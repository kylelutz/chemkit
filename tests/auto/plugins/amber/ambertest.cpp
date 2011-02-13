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
#include <chemkit/forcefield.h>
#include <chemkit/chemicalfile.h>
#include <chemkit/forcefieldatom.h>

const QString dataPath = "../../../data/";

class AmberTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void adenosine();
        void serine();
        void water();
};

void AmberTest::initTestCase()
{
    QVERIFY(chemkit::ForceField::forceFields().contains("amber"));
}

void AmberTest::adenosine()
{
    chemkit::ChemicalFile file(dataPath + "adenosine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), QString("C10H13N5O4"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->addMolecule(molecule);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    QCOMPARE(forceField->atomCount(), 32);
    QList<chemkit::ForceFieldAtom *> atoms = forceField->atoms();
    QCOMPARE(atoms[0]->type(), QString("CT"));
    QCOMPARE(atoms[1]->type(), QString("OS"));
    QCOMPARE(atoms[2]->type(), QString("CT"));
    QCOMPARE(atoms[3]->type(), QString("CT"));
    QCOMPARE(atoms[4]->type(), QString("OH"));
    QCOMPARE(atoms[5]->type(), QString("CT"));
    QCOMPARE(atoms[6]->type(), QString("OH"));
    QCOMPARE(atoms[7]->type(), QString("CT"));
    QCOMPARE(atoms[8]->type(), QString("OH"));
    QCOMPARE(atoms[9]->type(), QString("N*"));
    QCOMPARE(atoms[10]->type(), QString("CK"));
    QCOMPARE(atoms[11]->type(), QString("NB"));
    QCOMPARE(atoms[12]->type(), QString("CB"));
    QCOMPARE(atoms[13]->type(), QString("CB"));
    QCOMPARE(atoms[14]->type(), QString("NC"));
    QCOMPARE(atoms[15]->type(), QString("CQ"));
    QCOMPARE(atoms[16]->type(), QString("NC"));
    QCOMPARE(atoms[17]->type(), QString("CA"));
    QCOMPARE(atoms[18]->type(), QString("N2"));
    QCOMPARE(atoms[19]->type(), QString("H2"));
    QCOMPARE(atoms[20]->type(), QString("H1"));
    QCOMPARE(atoms[21]->type(), QString("H1"));
    QCOMPARE(atoms[22]->type(), QString("HO"));
    QCOMPARE(atoms[23]->type(), QString("H1"));
    QCOMPARE(atoms[24]->type(), QString("HO"));
    QCOMPARE(atoms[25]->type(), QString("H1"));
    QCOMPARE(atoms[26]->type(), QString("H1"));
    QCOMPARE(atoms[27]->type(), QString("HO"));
    QCOMPARE(atoms[28]->type(), QString("H5"));
    QCOMPARE(atoms[29]->type(), QString("H5"));
    QCOMPARE(atoms[30]->type(), QString("H"));
    QCOMPARE(atoms[31]->type(), QString("H"));

    QCOMPARE(forceField->calculationCount(), 585);
    QCOMPARE(qRound(forceField->energy()), 1460);

    delete forceField;
}

void AmberTest::serine()
{
    chemkit::ChemicalFile file(dataPath + "serine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), QString("C3H7NO3"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->addMolecule(molecule);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    QCOMPARE(forceField->atomCount(), 14);
    QList<chemkit::ForceFieldAtom *> atoms = forceField->atoms();
    QCOMPARE(atoms[0]->type(), QString("N3"));
    QCOMPARE(atoms[1]->type(), QString("CT"));
    QCOMPARE(atoms[2]->type(), QString("H"));
    QCOMPARE(atoms[3]->type(), QString("HP"));
    QCOMPARE(atoms[4]->type(), QString("C"));
    QCOMPARE(atoms[5]->type(), QString("CT"));
    QCOMPARE(atoms[6]->type(), QString("H1"));
    QCOMPARE(atoms[7]->type(), QString("H1"));
    QCOMPARE(atoms[8]->type(), QString("OH"));
    QCOMPARE(atoms[9]->type(), QString("HO"));
    QCOMPARE(atoms[10]->type(), QString("O2"));
    QCOMPARE(atoms[11]->type(), QString("O2"));
    QCOMPARE(atoms[12]->type(), QString("H"));
    QCOMPARE(atoms[13]->type(), QString("H"));

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
    QCOMPARE(forceField->atoms()[0]->type(), QString("OW"));
    QCOMPARE(forceField->atoms()[1]->type(), QString("HW"));
    QCOMPARE(forceField->atoms()[2]->type(), QString("HW"));

    QCOMPARE(forceField->calculationCount(), 3);

    QCOMPARE(qRound(forceField->energy()), 21085);

    delete forceField;
}

QTEST_APPLESS_MAIN(AmberTest)
#include "ambertest.moc"
