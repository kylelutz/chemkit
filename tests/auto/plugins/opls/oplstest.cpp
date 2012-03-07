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

// The OplsTest class validates the OPLS force field implementation
// in the opls plugin. Energies were checked against those
// calculated by TINKER.

#include "oplstest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculardescriptor.h>

const std::string dataPath = "../../../data/";

void OplsTest::initTestCase()
{
    // verify that the opls plugin registered itself correctly
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "opls") == 1);
    QVERIFY(boost::count(chemkit::ForceField::forceFields(), "opls") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "opls-energy") == 1);
}

void OplsTest::energy_data()
{
    QTest::addColumn<QString>("fileNameString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("energy");

    QTest::newRow("water") << "water.mol" << "H2O" << 1.8698;
    QTest::newRow("methanol") << "methanol.sdf" << "CH4O" << 5.6693;
    QTest::newRow("ethanol") << "ethanol.cml" << "C2H6O" << 1.2309;

}

void OplsTest::energy()
{
    QFETCH(QString, fileNameString);
    QFETCH(QString, formulaString);
    QFETCH(double, energy);

    QByteArray fileName = fileNameString.toAscii();
    QByteArray formula = formulaString.toAscii();

    boost::shared_ptr<chemkit::Molecule> molecule =
        chemkit::MoleculeFile::quickRead(dataPath + fileName.constData());
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula().c_str(), formula.constData());

    chemkit::ForceField *opls = chemkit::ForceField::create("opls");
    QVERIFY(opls != 0);

    opls->setMolecule(molecule.get());
    bool setup = opls->setup();
    QVERIFY(setup == true);

    QCOMPARE(qRound(opls->energy()), qRound(energy));

    // check opls energy descriptor
    QCOMPARE(qRound(molecule->descriptor("opls-energy").toDouble()), qRound(energy));

    delete opls;
}

QTEST_APPLESS_MAIN(OplsTest)
