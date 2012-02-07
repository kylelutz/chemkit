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

#include "ambertest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>
#include <chemkit/forcefieldatom.h>
#include <chemkit/moleculardescriptor.h>

const std::string dataPath = "../../../data/";

void AmberTest::initTestCase()
{
    // verify that the amber plugin registered itself correctly
    QVERIFY(boost::count(chemkit::ForceField::forceFields(), "amber") == 1);
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "amber") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "amber-energy") == 1);
}

void AmberTest::adenosine()
{
    chemkit::MoleculeFile file(dataPath + "adenosine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("C10H13N5O4"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->setMolecule(molecule.get());
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

    // check amber energy descriptor
    QCOMPARE(qRound(molecule->descriptor("amber-energy").toDouble()), 1460);

    delete forceField;
}

void AmberTest::serine()
{
    chemkit::MoleculeFile file(dataPath + "serine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("C3H7NO3"));

    chemkit::ForceField *forceField = chemkit::ForceField::create("amber");
    QVERIFY(forceField != 0);

    forceField->setMolecule(molecule.get());
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

    // check amber energy descriptor
    QCOMPARE(qRound(molecule->descriptor("amber-energy").toDouble()), 322);

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

    forceField->setMolecule(&water);
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
