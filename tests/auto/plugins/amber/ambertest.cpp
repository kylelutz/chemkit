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
#include <chemkit/topology.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculardescriptor.h>

#ifdef CHEMKIT_WITH_MD_IO
#include <chemkit/trajectoryfileformat.h>
#endif

const std::string dataPath = "../../../data/";

void AmberTest::initTestCase()
{
    // verify that the amber plugin registered itself correctly
    QVERIFY(boost::count(chemkit::ForceField::forceFields(), "amber") == 1);
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "amber") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "amber-energy") == 1);

    #ifdef CHEMKIT_WITH_MD_IO
    QVERIFY(boost::count(chemkit::TrajectoryFileFormat::formats(), "mdcrd") == 1);
    QVERIFY(boost::count(chemkit::TrajectoryFileFormat::formats(), "trj") == 1);
    #endif
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

    forceField->setTopologyFromMolecule(molecule.get());
    forceField->setup();
    QVERIFY(forceField->isSetup());

    boost::shared_ptr<chemkit::Topology> topology = forceField->topology();
    QCOMPARE(topology->size(), size_t(32));
    QCOMPARE(topology->type(0), std::string("CT"));
    QCOMPARE(topology->type(1), std::string("OS"));
    QCOMPARE(topology->type(2), std::string("CT"));
    QCOMPARE(topology->type(3), std::string("CT"));
    QCOMPARE(topology->type(4), std::string("OH"));
    QCOMPARE(topology->type(5), std::string("CT"));
    QCOMPARE(topology->type(6), std::string("OH"));
    QCOMPARE(topology->type(7), std::string("CT"));
    QCOMPARE(topology->type(8), std::string("OH"));
    QCOMPARE(topology->type(9), std::string("N*"));
    QCOMPARE(topology->type(10), std::string("CK"));
    QCOMPARE(topology->type(11), std::string("NB"));
    QCOMPARE(topology->type(12), std::string("CB"));
    QCOMPARE(topology->type(13), std::string("CB"));
    QCOMPARE(topology->type(14), std::string("NC"));
    QCOMPARE(topology->type(15), std::string("CQ"));
    QCOMPARE(topology->type(16), std::string("NC"));
    QCOMPARE(topology->type(17), std::string("CA"));
    QCOMPARE(topology->type(18), std::string("N2"));
    QCOMPARE(topology->type(19), std::string("H2"));
    QCOMPARE(topology->type(20), std::string("H1"));
    QCOMPARE(topology->type(21), std::string("H1"));
    QCOMPARE(topology->type(22), std::string("HO"));
    QCOMPARE(topology->type(23), std::string("H1"));
    QCOMPARE(topology->type(24), std::string("HO"));
    QCOMPARE(topology->type(25), std::string("H1"));
    QCOMPARE(topology->type(26), std::string("H1"));
    QCOMPARE(topology->type(27), std::string("HO"));
    QCOMPARE(topology->type(28), std::string("H5"));
    QCOMPARE(topology->type(29), std::string("H5"));
    QCOMPARE(topology->type(30), std::string("H"));
    QCOMPARE(topology->type(31), std::string("H"));

    QCOMPARE(forceField->calculationCount(), size_t(585));
    QCOMPARE(qRound(forceField->energy(molecule->coordinates())), 1460);

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

    forceField->setTopologyFromMolecule(molecule.get());
    forceField->setup();
    QVERIFY(forceField->isSetup());

    boost::shared_ptr<chemkit::Topology> topology = forceField->topology();
    QCOMPARE(topology->size(), size_t(14));
    QCOMPARE(topology->type(0), std::string("N3"));
    QCOMPARE(topology->type(1), std::string("CT"));
    QCOMPARE(topology->type(2), std::string("H"));
    QCOMPARE(topology->type(3), std::string("HP"));
    QCOMPARE(topology->type(4), std::string("C"));
    QCOMPARE(topology->type(5), std::string("CT"));
    QCOMPARE(topology->type(6), std::string("H1"));
    QCOMPARE(topology->type(7), std::string("H1"));
    QCOMPARE(topology->type(8), std::string("OH"));
    QCOMPARE(topology->type(9), std::string("HO"));
    QCOMPARE(topology->type(10), std::string("O2"));
    QCOMPARE(topology->type(11), std::string("O2"));
    QCOMPARE(topology->type(12), std::string("H"));
    QCOMPARE(topology->type(13), std::string("H"));

    QCOMPARE(forceField->calculationCount(), size_t(118));
    QCOMPARE(qRound(forceField->energy(molecule->coordinates())), 322);

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

    forceField->setTopologyFromMolecule(&water);
    forceField->setup();
    QVERIFY(forceField->isSetup());

    boost::shared_ptr<chemkit::Topology> topology = forceField->topology();
    QCOMPARE(topology->size(), size_t(3));
    QCOMPARE(topology->type(0), std::string("OW"));
    QCOMPARE(topology->type(1), std::string("HW"));
    QCOMPARE(topology->type(2), std::string("HW"));

    QCOMPARE(forceField->calculationCount(), size_t(3));

    QCOMPARE(qRound(forceField->energy(water.coordinates())), 21085);

    delete forceField;
}

QTEST_APPLESS_MAIN(AmberTest)
