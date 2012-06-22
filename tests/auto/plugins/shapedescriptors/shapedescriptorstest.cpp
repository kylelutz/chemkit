/******************************************************************************
**
** Copyright (C) 2012 Kitware, Inc.
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

#include "shapedescriptorstest.h"

#include <boost/shared_ptr.hpp>
#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculardescriptor.h>

const std::string dataPath = "../../../data/";

void ShapeDescriptorsTest::initTestCase()
{
    // verify that the shapedescriptors plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "gravitational-index") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "bonded-gravitational-index") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "radius-of-gyration") == 1);
}

void ShapeDescriptorsTest::test_data()
{
    QTest::addColumn<QString>("fileNameString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("gravitationalIndex");
    QTest::addColumn<double>("bondedGravitationalIndex");
    QTest::addColumn<double>("radiusOfGyration");

    QTest::newRow("ethanol") << "ethanol.cml" << "C2H6O" << 295.156 << 227.184 << 1.22665;
    QTest::newRow("glucose") << "glucose.cml" << "C6H12O6" << 2399.91 << 1126.53 << 2.41039;
    QTest::newRow("guanine") << "guanine.cml" << "C5H5N5O" << 2117.31 << 1159.68 << 2.19339;
    QTest::newRow("serine") << "serine.mol" << "C3H7NO3" << 1102.91 << 623.636 << 2.6157;
    QTest::newRow("uridine") << "uridine.mol2" << "C9H13N2O9P" << 5417.02 << 2499.18 << 3.47502;
}

void ShapeDescriptorsTest::test()
{
    QFETCH(QString, fileNameString);
    QFETCH(QString, formulaString);
    QFETCH(double, gravitationalIndex);
    QFETCH(double, bondedGravitationalIndex);
    QFETCH(double, radiusOfGyration);

    QByteArray fileName = fileNameString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::MoleculeFile file(dataPath + fileName.constData());
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula().c_str(), formula.constData());
    QCOMPARE(qRound(molecule->descriptor("gravitational-index").toDouble()), qRound(gravitationalIndex));
    QCOMPARE(qRound(molecule->descriptor("bonded-gravitational-index").toDouble()), qRound(bondedGravitationalIndex));
    QCOMPARE(qRound(molecule->descriptor("radius-of-gyration").toDouble() * 100), qRound(radiusOfGyration * 100));
}

QTEST_APPLESS_MAIN(ShapeDescriptorsTest)
