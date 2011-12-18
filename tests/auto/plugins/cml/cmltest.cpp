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

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/coordinateset.h>
#include <chemkit/diagramcoordinates.h>
#include <chemkit/moleculefileformat.h>
#include <chemkit/cartesiancoordinates.h>

const std::string dataPath = "../../../data/";

void CmlTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::MoleculeFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "cml") != formats.end());
}

void CmlTest::read_data()
{
    QTest::addColumn<QString>("fileName");
    QTest::addColumn<QString>("formula");

    QTest::newRow("buckminsterfullerene") << "buckminsterfullerene.cml" << "C60";
    QTest::newRow("ethanol") << "ethanol.cml" << "C2H6O";
    QTest::newRow("gaunine") << "guanine.cml" << "C5H5N5O";
    QTest::newRow("paracetamol") << "paracetamol.cml" << "C8H9NO2";
}

void CmlTest::read()
{
    QFETCH(QString, fileName);
    QFETCH(QString, formula);

    chemkit::MoleculeFile file(dataPath + fileName.toStdString());
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), formula.toStdString());
}

void CmlTest::glucose()
{
    chemkit::MoleculeFile file(dataPath + "glucose.cml");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C6H12O6"));

    QCOMPARE(molecule->coordinateSetCount(), size_t(2));
    QVERIFY(molecule->coordinateSet(0)->type() == chemkit::CoordinateSet::Cartesian);
    QVERIFY(molecule->coordinateSet(1)->type() == chemkit::CoordinateSet::Diagram);

    chemkit::CartesianCoordinates *cartesianCoordinates = molecule->coordinateSet(0)->cartesianCoordinates();
    QVERIFY(cartesianCoordinates != 0);
    QCOMPARE(cartesianCoordinates->size(), 24);

    chemkit::DiagramCoordinates *diagramCoordinates = molecule->coordinateSet(1)->diagramCoordinates();
    QVERIFY(diagramCoordinates != 0);
    QCOMPARE(diagramCoordinates->size(), size_t(24));
}

QTEST_APPLESS_MAIN(CmlTest)
