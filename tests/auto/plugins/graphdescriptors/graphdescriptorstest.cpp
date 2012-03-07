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

#include "graphdescriptorstest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void GraphDescriptorsTest::initTestCase()
{
    // verify that the graphdescriptors plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "graph-density") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "graph-diameter") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "graph-order") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "graph-radius") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "graph-size") == 1);
}

void GraphDescriptorsTest::test_data()
{
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("graphDensity");
    QTest::addColumn<int>("graphDiameter");
    QTest::addColumn<int>("graphOrder");
    QTest::addColumn<int>("graphRadius");
    QTest::addColumn<int>("graphSize");

    QTest::newRow("ethane") << "CC" << "C2H6" << 0.25 << 3 << 8 << 2 << 7;
    QTest::newRow("ethanol") << "CCO" << "C2H6O" << 0.222222 << 4 << 9 << 2 << 8;
    QTest::newRow("butane") << "CCC" << "C3H8" << 0.181818 << 4 << 11 << 2 << 10;
    QTest::newRow("cyclohexane") << "C1CCCCC1" << "C6H12" << 0.117647 << 5 << 18 << 4 << 18;
    QTest::newRow("toluene") << "Cc1ccccc1" << "C7H8" << 0.142857 << 6 << 15 << 4 << 15;
    QTest::newRow("naphthalene") << "c1ccc2ccccc2c1" << "C10H8" << 0.124183 << 7 << 18 << 4 << 19;
    QTest::newRow("biotin") << "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12" << "C10H16N2O3S" << 0.0665323 << 11 << 32 << 6 << 33;
    QTest::newRow("adenosine") << "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O" << "C10H13N5O4" << 0.0685484 << 11 << 32 << 6 << 34;
}

void GraphDescriptorsTest::test()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);
    QFETCH(double, graphDensity);
    QFETCH(int, graphDiameter);
    QFETCH(int, graphOrder);
    QFETCH(int, graphRadius);
    QFETCH(int, graphSize);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());
    QCOMPARE(qRound(molecule.descriptor("graph-density").toDouble() * 100.0), qRound(graphDensity * 100.0));
    QCOMPARE(molecule.descriptor("graph-diameter").toInt(), graphDiameter);
    QCOMPARE(molecule.descriptor("graph-order").toInt(), graphOrder);
    QCOMPARE(molecule.descriptor("graph-radius").toInt(), graphRadius);
    QCOMPARE(molecule.descriptor("graph-size").toInt(), graphSize);
}

QTEST_APPLESS_MAIN(GraphDescriptorsTest)
