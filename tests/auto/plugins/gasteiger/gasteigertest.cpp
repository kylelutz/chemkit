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
#include <chemkit/partialchargepredictor.h>

class GasteigerTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void name();
        void methane();
        void fluoromethane();
        void ethane();
        void fluoroethane();
};

void GasteigerTest::initTestCase()
{
    QVERIFY(chemkit::PartialChargePredictor::predictors().contains("gasteiger"));
}

void GasteigerTest::name()
{
    chemkit::PartialChargePredictor *predictor = chemkit::PartialChargePredictor::create("gasteiger");
    QVERIFY(predictor != 0);
    QCOMPARE(predictor->name(), QString("gasteiger"));

    delete predictor;
}

void GasteigerTest::methane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Atom *H3 = molecule.addAtom("H");
    chemkit::Atom *H4 = molecule.addAtom("H");
    chemkit::Atom *H5 = molecule.addAtom("H");
    molecule.addBond(C1, H2);
    molecule.addBond(C1, H3);
    molecule.addBond(C1, H4);
    molecule.addBond(C1, H5);
    QCOMPARE(molecule.formula(), std::string("CH4"));

    chemkit::PartialChargePredictor *predictor = chemkit::PartialChargePredictor::create("gasteiger");
    QVERIFY(predictor != 0);

    predictor->setMolecule(&molecule);
    QCOMPARE(qRound(predictor->partialCharge(C1) * 1e3), -78);

    delete predictor;
}

void GasteigerTest::fluoromethane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *F2 = molecule.addAtom("F");
    chemkit::Atom *H3 = molecule.addAtom("H");
    chemkit::Atom *H4 = molecule.addAtom("H");
    chemkit::Atom *H5 = molecule.addAtom("H");
    molecule.addBond(C1, F2);
    molecule.addBond(C1, H3);
    molecule.addBond(C1, H4);
    molecule.addBond(C1, H5);
    QCOMPARE(molecule.formula(), std::string("CH3F"));

    chemkit::PartialChargePredictor *predictor = chemkit::PartialChargePredictor::create("gasteiger");
    QVERIFY(predictor != 0);

    predictor->setMolecule(&molecule);
    QCOMPARE(qRound(predictor->partialCharge(C1) * 1e3), 79);
    QCOMPARE(qRound(predictor->partialCharge(F2) * 1e3), -253);
    QCOMPARE(qRound(predictor->partialCharge(H3) * 1e3), 58);
    QCOMPARE(qRound(predictor->partialCharge(H4) * 1e3), 58);
    QCOMPARE(qRound(predictor->partialCharge(H5) * 1e3), 58);

    delete predictor;
}

void GasteigerTest::ethane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Atom *H3 = molecule.addAtom("H");
    chemkit::Atom *H4 = molecule.addAtom("H");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *H6 = molecule.addAtom("H");
    chemkit::Atom *H7 = molecule.addAtom("H");
    chemkit::Atom *H8 = molecule.addAtom("H");
    molecule.addBond(C1, H2);
    molecule.addBond(C1, H3);
    molecule.addBond(C1, H4);
    molecule.addBond(C1, C5);
    molecule.addBond(C5, H6);
    molecule.addBond(C5, H7);
    molecule.addBond(C5, H8);
    QCOMPARE(molecule.formula(), std::string("C2H6"));

    chemkit::PartialChargePredictor *predictor = chemkit::PartialChargePredictor::create("gasteiger");
    QVERIFY(predictor != 0);

    predictor->setMolecule(&molecule);
    QCOMPARE(qRound(predictor->partialCharge(C1) * 1e3), -68);
    QCOMPARE(qRound(predictor->partialCharge(C5) * 1e3), -68);

    delete predictor;
}

void GasteigerTest::fluoroethane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Atom *H3 = molecule.addAtom("H");
    chemkit::Atom *H4 = molecule.addAtom("H");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *F6 = molecule.addAtom("F");
    chemkit::Atom *H7 = molecule.addAtom("H");
    chemkit::Atom *H8 = molecule.addAtom("H");
    molecule.addBond(C1, H2);
    molecule.addBond(C1, H3);
    molecule.addBond(C1, H4);
    molecule.addBond(C1, C5);
    molecule.addBond(C5, F6);
    molecule.addBond(C5, H7);
    molecule.addBond(C5, H8);
    QCOMPARE(molecule.formula(), std::string("C2H5F"));

    chemkit::PartialChargePredictor *predictor = chemkit::PartialChargePredictor::create("gasteiger");
    QVERIFY(predictor != 0);

    predictor->setMolecule(&molecule);
    QCOMPARE(qRound(predictor->partialCharge(C1) * 1e3), -37);
    QCOMPARE(qRound(predictor->partialCharge(C5) * 1e3), 87);

    delete predictor;
}

QTEST_APPLESS_MAIN(GasteigerTest)
#include "gasteigertest.moc"
