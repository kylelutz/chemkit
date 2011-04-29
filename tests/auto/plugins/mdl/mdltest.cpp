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

#include "mdltest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculefileformat.h>

const std::string dataPath = "../../../data/";

void MdlTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::MoleculeFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "mol") != formats.end());
    QVERIFY(std::find(formats.begin(), formats.end(), "mdl") != formats.end());
    QVERIFY(std::find(formats.begin(), formats.end(), "sdf") != formats.end());
    QVERIFY(std::find(formats.begin(), formats.end(), "sd") != formats.end());
}

void MdlTest::read_methanol()
{
    chemkit::MoleculeFile file(dataPath + "methanol.sdf");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // check molecule
    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("CH4O"));

    // check data
    QCOMPARE(molecule->name(), std::string("887"));
    QCOMPARE(molecule->data("PUBCHEM_COMPOUND_CID").toString(), QString("887"));
    QCOMPARE(molecule->data("PUBCHEM_HEAVY_ATOM_COUNT").toInt(), 2);
}

void MdlTest::read_guanine()
{
    chemkit::MoleculeFile file(dataPath + "guanine.mol");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QCOMPARE(ok, true);

    // check format
    QVERIFY(file.format() != 0);
    QCOMPARE(file.formatName(), std::string("mol"));

    // check molecule
    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *guanine = file.molecule();
    QCOMPARE(guanine->formula(), std::string("C5H5N5O"));
    QCOMPARE(guanine->name(), std::string("Guanine"));
    QCOMPARE(guanine->atomCount(), 16);
    QCOMPARE(guanine->bondCount(), 17);
}

void MdlTest::read_benzenes()
{
    chemkit::MoleculeFile file(dataPath + "pubchem_416_benzenes.sdf");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QCOMPARE(ok, true);

    // check format
    QVERIFY(file.format() != 0);
    QCOMPARE(file.formatName(), std::string("sdf"));

    // check molecules
    QCOMPARE(file.moleculeCount(), 416);

    // check molecule data
    foreach(const chemkit::Molecule *molecule, file.molecules()){
        QCOMPARE(molecule->name(), molecule->data("PUBCHEM_COMPOUND_CID").toString().toStdString());
    }
}

QTEST_APPLESS_MAIN(MdlTest)
