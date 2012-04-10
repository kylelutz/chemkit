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

#include "gasteigertest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/partialchargemodel.h>

void GasteigerTest::initTestCase()
{
    // verify that the gasteiger plugin registered itself correctly
    QVERIFY(boost::count(chemkit::PartialChargeModel::models(), "gasteiger") == 1);
}

void GasteigerTest::name()
{
    chemkit::PartialChargeModel *model = chemkit::PartialChargeModel::create("gasteiger");
    QVERIFY(model != 0);
    QCOMPARE(model->name(), std::string("gasteiger"));

    delete model;
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

    chemkit::PartialChargeModel *model = chemkit::PartialChargeModel::create("gasteiger");
    QVERIFY(model != 0);

    model->setMolecule(&molecule);
    QCOMPARE(qRound(model->partialCharge(C1) * 1e3), -78);

    delete model;
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

    chemkit::PartialChargeModel *model = chemkit::PartialChargeModel::create("gasteiger");
    QVERIFY(model != 0);

    model->setMolecule(&molecule);
    QCOMPARE(qRound(model->partialCharge(C1) * 1e3), 79);
    QCOMPARE(qRound(model->partialCharge(F2) * 1e3), -253);
    QCOMPARE(qRound(model->partialCharge(H3) * 1e3), 58);
    QCOMPARE(qRound(model->partialCharge(H4) * 1e3), 58);
    QCOMPARE(qRound(model->partialCharge(H5) * 1e3), 58);

    delete model;
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

    chemkit::PartialChargeModel *model = chemkit::PartialChargeModel::create("gasteiger");
    QVERIFY(model != 0);

    model->setMolecule(&molecule);
    QCOMPARE(qRound(model->partialCharge(C1) * 1e3), -68);
    QCOMPARE(qRound(model->partialCharge(C5) * 1e3), -68);

    delete model;
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

    chemkit::PartialChargeModel *model = chemkit::PartialChargeModel::create("gasteiger");
    QVERIFY(model != 0);

    model->setMolecule(&molecule);
    QCOMPARE(qRound(model->partialCharge(C1) * 1e3), -37);
    QCOMPARE(qRound(model->partialCharge(C5) * 1e3), 87);

    delete model;
}

QTEST_APPLESS_MAIN(GasteigerTest)
