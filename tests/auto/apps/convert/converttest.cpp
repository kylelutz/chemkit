/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include <QtTest>

#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>

const QString testDataPath = "../../../data/";
const QString convertApplication = "../../../../bin/chemkit-convert";

// The ConvertTest class tests the chemkit-convert tool. Each test
// case converts a file to a different file format. The molecules in
// the converted file are then verified against the original file to
// ensure the conversion was correct. Each molecule is verified using
// the Molecule::name() and Molecule::equals() methods.
class ConvertTest : public QObject
{
    Q_OBJECT

    private slots:
        void convertEthanol();
        void convertBenzenes();
};

void ConvertTest::convertEthanol()
{
    // read input file
    chemkit::ChemicalFile inputFile(testDataPath + "ethanol.cml");
    bool ok = inputFile.read();
    if(!ok){
        qDebug() << inputFile.errorString();
    }
    QCOMPARE(ok, true);
    QCOMPARE(inputFile.moleculeCount(), 1);
    chemkit::Molecule *inputEthanol = inputFile.molecule();
    QCOMPARE(inputEthanol->formula(), QString("C2H6O"));

    // setup output file
    QTemporaryFile output("XXXXXX.mol");
    output.open();

    // setup arguments
    QStringList arguments;
    arguments.append(inputFile.fileName());
    arguments.append(output.fileName());

    // run chemkit-convert
    QProcess process;
    process.start(convertApplication, arguments);
    process.waitForFinished();
    process.close();

    // read and verify the output file
    chemkit::ChemicalFile outputFile;
    ok = outputFile.read(output.fileName());
    if(!ok){
        qDebug() << outputFile.errorString();
    }
    QCOMPARE(ok, true);
    QCOMPARE(outputFile.moleculeCount(), 1);

    // verify the output molecule
    chemkit::Molecule *outputEthanol = outputFile.molecule();
    QCOMPARE(outputEthanol->name(), inputEthanol->name());
    QCOMPARE(outputEthanol->equals(inputEthanol), true);
}

void ConvertTest::convertBenzenes()
{
    // read input file
    chemkit::ChemicalFile inputFile(testDataPath + "pubchem_416_benzenes.sdf");
    bool ok = inputFile.read();
    if(!ok){
        qDebug() << inputFile.errorString();
    }
    QCOMPARE(ok, true);
    QCOMPARE(inputFile.moleculeCount(), 416);

    // setup output file
    QTemporaryFile tempFile("XXXXXX.mol2");
    tempFile.open();

    // setup arguments
    QStringList arguments;
    arguments.append(inputFile.fileName());
    arguments.append(tempFile.fileName());

    // run chemkit-convert
    QProcess process;
    process.start(convertApplication, arguments);
    process.waitForFinished();
    process.close();

    // read and verify the output file
    chemkit::ChemicalFile outputFile;
    ok = outputFile.read(tempFile.fileName());
    if(!ok){
        qDebug() << outputFile.errorString();
    }
    QCOMPARE(ok, true);
    QCOMPARE(outputFile.moleculeCount(), 416);

    // verify the output molecules
    for(int i = 0; i < 416; i++){
        chemkit::Molecule *inputMolecule = inputFile.molecule(i);
        chemkit::Molecule *outputMolecule = outputFile.molecule(i);

        QCOMPARE(inputMolecule->name(), outputMolecule->name());
        QCOMPARE(inputMolecule->equals(outputMolecule, chemkit::Molecule::CompareAromaticity | chemkit::Molecule::CompareHydrogens), true);
    }
}

QTEST_APPLESS_MAIN(ConvertTest)
#include "converttest.moc"
