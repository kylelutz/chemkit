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

#include "greptest.h"
 
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const QString grepApplication = "../../../../bin/chemkit-grep";
const QString testDataPath = "../../../data/";

void GrepTest::ironComposition()
{
    QProcess process;

    QStringList arguments;
    arguments.append("--composition");
    arguments.append("InChI=1/Fe/");
    arguments.append(testDataPath + "MMFF94_hypervalent.mol2");

    QTemporaryFile output;
    output.open();

    process.setStandardOutputFile(output.fileName());
    process.start(grepApplication, arguments);
    process.waitForFinished();
    process.close();

    chemkit::MoleculeFile file;
    QByteArray outputFileName = output.fileName().toAscii();
    bool ok = file.read(outputFileName.constData(), "mol2");
    if(!ok){
        qDebug() << file.errorString().c_str();
    }
    QVERIFY(ok);
    QCOMPARE(file.moleculeCount(), size_t(2));
    QCOMPARE(file.molecule(0)->name(), std::string("FE2PW3"));
    QCOMPARE(file.molecule(0)->formula(), std::string("FeH6O3"));
    QCOMPARE(file.molecule(1)->name(), std::string("FE3PW3"));
    QCOMPARE(file.molecule(1)->formula(), std::string("FeH6O3"));
}

void GrepTest::benzoicAcid()
{
    QProcess process;

    QStringList arguments;
    arguments.append("--names-only");
    arguments.append("InChI=1/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)/f/h8H");
    arguments.append(testDataPath + "pubchem_416_benzenes.sdf");

    process.start(grepApplication, arguments);
    process.waitForFinished();

    QString output = process.readAllStandardOutput();
    QStringList moleculeNames = output.split("\n", QString::SkipEmptyParts);
    QCOMPARE(moleculeNames.size(), 39);

    // trim whitespace from names
    for(int i = 0; i < moleculeNames.size(); i++){
        moleculeNames[i] = moleculeNames[i].trimmed();
    }

    QCOMPARE(moleculeNames[0], QString("2605"));
    QCOMPARE(moleculeNames[1], QString("2541"));
    QCOMPARE(moleculeNames[2], QString("2536"));
    QCOMPARE(moleculeNames[3], QString("2480"));
    QCOMPARE(moleculeNames[4], QString("2471"));
    QCOMPARE(moleculeNames[5], QString("2424"));
    QCOMPARE(moleculeNames[6], QString("2418"));
    QCOMPARE(moleculeNames[7], QString("2376"));
    QCOMPARE(moleculeNames[8], QString("2356"));
    QCOMPARE(moleculeNames[9], QString("2347"));
    QCOMPARE(moleculeNames[10], QString("2334"));
    QCOMPARE(moleculeNames[11], QString("2329"));
    QCOMPARE(moleculeNames[12], QString("2318"));
    QCOMPARE(moleculeNames[13], QString("2316"));
    QCOMPARE(moleculeNames[14], QString("2259"));
    QCOMPARE(moleculeNames[15], QString("2258"));
    QCOMPARE(moleculeNames[16], QString("2257"));
    QCOMPARE(moleculeNames[17], QString("2126"));
    QCOMPARE(moleculeNames[18], QString("1974"));
    QCOMPARE(moleculeNames[19], QString("1964"));
    QCOMPARE(moleculeNames[20], QString("1854"));
    QCOMPARE(moleculeNames[21], QString("1829"));
    QCOMPARE(moleculeNames[22], QString("1730"));
    QCOMPARE(moleculeNames[23], QString("1718"));
    QCOMPARE(moleculeNames[24], QString("1711"));
    QCOMPARE(moleculeNames[25], QString("1683"));
    QCOMPARE(moleculeNames[26], QString("1583"));
    QCOMPARE(moleculeNames[27], QString("1559"));
    QCOMPARE(moleculeNames[28], QString("1385"));
    QCOMPARE(moleculeNames[29], QString("1311"));
    QCOMPARE(moleculeNames[30], QString("1309"));
    QCOMPARE(moleculeNames[31], QString("1287"));
    QCOMPARE(moleculeNames[32], QString("1261"));
    QCOMPARE(moleculeNames[33], QString("1248"));
    QCOMPARE(moleculeNames[34], QString("894"));
    QCOMPARE(moleculeNames[35], QString("514"));
    QCOMPARE(moleculeNames[36], QString("287"));
    QCOMPARE(moleculeNames[37], QString("341"));
    QCOMPARE(moleculeNames[38], QString("170"));
}

QTEST_APPLESS_MAIN(GrepTest)
