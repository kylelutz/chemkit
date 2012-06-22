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

#include "chemjsontest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/moleculefileformat.h>

const std::string dataPath = "../../../data/";

void ChemJsonTest::initTestCase()
{
    // verify that the chemjson plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "cjson") == 1);
}

void ChemJsonTest::read()
{
    chemkit::MoleculeFile file(dataPath + "ethane.cjson");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.moleculeCount(), size_t(1));

    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), std::string("C2H6"));
    QCOMPARE(molecule->atomCount(), size_t(8));
    QCOMPARE(molecule->bondCount(), size_t(7));
    QCOMPARE(molecule->data("molecular weight").toInt(), 30);
    QCOMPARE(molecule->data("melting point").toInt(), -172);
    QCOMPARE(molecule->data("boiling point").toInt(), -88);
}

void ChemJsonTest::checkSanitizeStrings()
{
    // create helium molecule
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    molecule->addAtom("He");
    QCOMPARE(molecule->formula(), std::string("He"));

    // set the molecule's name to a string containing quotes
    molecule->setName("\"helium\"");

    // create file
    chemkit::MoleculeFile outputFile;
    outputFile.addMolecule(molecule);

    // write file
    std::stringstream stream;
    bool ok = outputFile.write(stream, "cjson");
    if(!ok)
        qDebug() << outputFile.errorString().c_str();
    QVERIFY(ok);

    // read file
    chemkit::MoleculeFile inputFile;
    stream.seekg(0);
    ok = inputFile.read(stream, "cjson");
    if(!ok)
        qDebug() << inputFile.errorString().c_str();
    QVERIFY(ok);

    // check the molecule
    molecule = inputFile.molecule();
    QCOMPARE(molecule->formula(), std::string("He"));
    QCOMPARE(molecule->name(), std::string("helium"));
}

QTEST_APPLESS_MAIN(ChemJsonTest)
