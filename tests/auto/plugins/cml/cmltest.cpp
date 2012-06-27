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

#include "cmltest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/coordinateset.h>
#include <chemkit/diagramcoordinates.h>
#include <chemkit/moleculefileformat.h>
#include <chemkit/cartesiancoordinates.h>

const std::string dataPath = "../../../data/";

void CmlTest::initTestCase()
{
    // verify that the cml plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "cml") == 1);
}

void CmlTest::read_data()
{
    QTest::addColumn<QString>("fileNameString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<int>("atomCount");
    QTest::addColumn<int>("bondCount");

    QTest::newRow("buckminsterfullerene") << "buckminsterfullerene.cml" << "C60" << 60 << 90;
    QTest::newRow("ethanol") << "ethanol.cml" << "C2H6O" << 9 << 8;
    QTest::newRow("gaunine") << "guanine.cml" << "C5H5N5O" << 16 << 17;
    QTest::newRow("paracetamol") << "paracetamol.cml" << "C8H9NO2" << 20 << 20;
}

void CmlTest::read()
{
    QFETCH(QString, fileNameString);
    QFETCH(QString, formulaString);
    QFETCH(int, atomCount);
    QFETCH(int, bondCount);

    QByteArray fileName = fileNameString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::MoleculeFile file(dataPath + fileName.constData());
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(1));
    const boost::shared_ptr<chemkit::Molecule> &molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula().c_str(), formula.constData());
    QCOMPARE(molecule->atomCount(), size_t(atomCount));
    QCOMPARE(molecule->bondCount(), size_t(bondCount));
    QCOMPARE(molecule->coordinateSetCount(), size_t(1));
    QVERIFY(molecule->coordinateSet(0)->type() == chemkit::CoordinateSet::Cartesian);
}

void CmlTest::glucose()
{
    chemkit::MoleculeFile file(dataPath + "glucose.cml");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C6H12O6"));

    QCOMPARE(molecule->coordinateSetCount(), size_t(2));
    QVERIFY(molecule->coordinateSet(0)->type() == chemkit::CoordinateSet::Cartesian);
    QVERIFY(molecule->coordinateSet(1)->type() == chemkit::CoordinateSet::Diagram);

    chemkit::CartesianCoordinates *cartesianCoordinates = molecule->coordinateSet(0)->cartesianCoordinates();
    QVERIFY(cartesianCoordinates != 0);
    QCOMPARE(cartesianCoordinates->size(), size_t(24));

    chemkit::DiagramCoordinates *diagramCoordinates = molecule->coordinateSet(1)->diagramCoordinates();
    QVERIFY(diagramCoordinates != 0);
    QCOMPARE(diagramCoordinates->size(), size_t(24));
}

void CmlTest::ethanol()
{
    chemkit::MoleculeFile file(dataPath + "ethanol.cml");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C2H6O"));

    // check molecule property data
    QCOMPARE(qRound(molecule->data("molecular weight").toReal()), 46);
    QCOMPARE(qRound(molecule->data("exact molecular weight").toReal()), 46);
    QCOMPARE(qRound(molecule->data("melting point").toReal()), -114);
    QCOMPARE(qRound(molecule->data("boiling point").toReal()), 78);
}

QTEST_APPLESS_MAIN(CmlTest)
