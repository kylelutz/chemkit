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

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculefileformat.h>

const std::string dataPath = "../../../data/";

void MdlTest::initTestCase()
{
    // verify that the mdl plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "mol") == 1);
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "mdl") == 1);
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "sdf") == 1);
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "sd") == 1);
}

void MdlTest::read_methanol()
{
    chemkit::MoleculeFile file(dataPath + "methanol.sdf");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // check molecule
    QCOMPARE(file.moleculeCount(), size_t(1));
    const boost::shared_ptr<chemkit::Molecule> &molecule = file.molecule();
    QCOMPARE(molecule->formula(), std::string("CH4O"));

    // check data
    QCOMPARE(molecule->name(), std::string("887"));
    QCOMPARE(molecule->data("PUBCHEM_COMPOUND_CID").toString(), std::string("887"));
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
    QCOMPARE(file.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> guanine = file.molecule();
    QCOMPARE(guanine->formula(), std::string("C5H5N5O"));
    QCOMPARE(guanine->name(), std::string("Guanine"));
    QCOMPARE(guanine->atomCount(), size_t(16));
    QCOMPARE(guanine->bondCount(), size_t(17));
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
    QCOMPARE(file.moleculeCount(), size_t(416));

    // check molecule data
    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, file.molecules()){
        QCOMPARE(molecule->name(), molecule->data("PUBCHEM_COMPOUND_CID").toString());
    }
}

void MdlTest::read_serine()
{
    // check that gz compression is supported
    const std::vector<std::string> &compressionFormats = chemkit::MoleculeFile::compressionFormats();
    if(std::find(compressionFormats.begin(),
                 compressionFormats.end(),
                 "gz") == compressionFormats.end()){
        QSKIP("Gzip compression not supported", SkipSingle);
    }

    chemkit::MoleculeFile file(dataPath + "serine.mol.gz");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // check format
    QVERIFY(file.format() != 0);
    QCOMPARE(file.formatName(), std::string("mol"));

    // check molecule
    const boost::shared_ptr<chemkit::Molecule> &molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C3H7NO3"));
}

QTEST_APPLESS_MAIN(MdlTest)
