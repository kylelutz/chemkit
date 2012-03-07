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

// The ConvertTest class tests the chemkit-convert tool. Each test
// case converts a file to a different file format. The molecules in
// the converted file are then verified against the original file to
// ensure the conversion was correct. Each molecule is verified using
// the Molecule::name() and Molecule::equals() methods.

#include "converttest.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const std::string testDataPath = "../../../data/";
const QString convertApplication = "../../../../bin/chemkit-convert";

void ConvertTest::convertEthanol()
{
    // read input file
    chemkit::MoleculeFile inputFile(testDataPath + "ethanol.cml");
    bool ok = inputFile.read();
    if(!ok){
        qDebug() << inputFile.errorString().c_str();
    }
    QCOMPARE(ok, true);
    QCOMPARE(inputFile.moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> inputEthanol = inputFile.molecule();
    QCOMPARE(inputEthanol->formula(), std::string("C2H6O"));

    // setup output file
    QTemporaryFile output("XXXXXX.mol");
    output.open();

    // setup arguments
    QStringList arguments;
    arguments.append(inputFile.fileName().c_str());
    arguments.append(output.fileName());

    // run chemkit-convert
    QProcess process;
    process.start(convertApplication, arguments);
    process.waitForFinished();
    process.close();

    // read and verify the output file
    QByteArray outputFileName = output.fileName().toAscii();
    chemkit::MoleculeFile outputFile;
    ok = outputFile.read(outputFileName.constData());
    if(!ok){
        qDebug() << outputFile.errorString().c_str();
    }
    QCOMPARE(ok, true);
    QCOMPARE(outputFile.moleculeCount(), size_t(1));

    // verify the output molecule
    boost::shared_ptr<chemkit::Molecule> outputEthanol = outputFile.molecule();
    QCOMPARE(outputEthanol->name(), inputEthanol->name());
}

void ConvertTest::convertBenzenes()
{
    // read input file
    chemkit::MoleculeFile inputFile(testDataPath + "pubchem_416_benzenes.sdf");
    bool ok = inputFile.read();
    if(!ok){
        qDebug() << inputFile.errorString().c_str();
    }
    QCOMPARE(ok, true);
    QCOMPARE(inputFile.moleculeCount(), size_t(416));

    // setup output file
    QTemporaryFile tempFile("XXXXXX.mol2");
    tempFile.open();

    // setup arguments
    QStringList arguments;
    arguments.append(inputFile.fileName().c_str());
    arguments.append(tempFile.fileName());

    // run chemkit-convert
    QProcess process;
    process.start(convertApplication, arguments);
    process.waitForFinished();
    process.close();

    // read and verify the output file
    chemkit::MoleculeFile outputFile;
    QByteArray outputFileName = tempFile.fileName().toAscii();
    ok = outputFile.read(outputFileName.constData());
    if(!ok){
        qDebug() << outputFile.errorString().c_str();
    }
    QCOMPARE(ok, true);
    QCOMPARE(outputFile.moleculeCount(), size_t(416));

    // verify the output molecules
    for(int i = 0; i < 416; i++){
        boost::shared_ptr<chemkit::Molecule> inputMolecule = inputFile.molecule(i);
        boost::shared_ptr<chemkit::Molecule> outputMolecule = outputFile.molecule(i);

        QCOMPARE(inputMolecule->name(), outputMolecule->name());
    }
}

QTEST_APPLESS_MAIN(ConvertTest)
