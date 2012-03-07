/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "pubchemtest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/fingerprint.h>

namespace {

// compare two fingerprints
void compareFingerprints(const chemkit::Bitset &actual, const boost::dynamic_bitset<unsigned char> &expected)
{
    for(size_t i = 0; i < actual.size(); i++){
        // skip section 3 - ring counts
        if(i >= 115 && i <= 262){
            continue;
        }

        // skip section 4 - simple atom nearest neighbors
        if(i >= 327 && i <= 415){
            continue;
        }

        // skip section 5 - detailed atom neighborhoods
        if(i >= 416 && i <= 459){
            continue;
        }

        // skip section 6 - simples SMARTS patterns
        if(i >= 460 && i <= 712){
            continue;
        }

        // skip section 7 - complex SMARTS patterns
        if(i >= 713 && i <= 880){
            continue;
        }

        if(actual[i] != expected[i]){
            qDebug() << "Value at bit" << i << "does not match.";
            qDebug() << "  Actual:" << actual[i];
            qDebug() << "  Expected:" << expected[i];
            QVERIFY(actual[i] == expected[i]);
        }
    }
}

// reverse bits in a byte
unsigned char reverse(unsigned char b)
{
   b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
   b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
   b = (b & 0xAA) >> 1 | (b & 0x55) << 1;

   return b;
}

} // end anonymous namespace

void PubChemTest::initTestCase()
{
    // verify that the pubchem plugin registered itself correctly
    QVERIFY(boost::count(chemkit::Fingerprint::fingerprints(), "pubchem") == 1);
}

void PubChemTest::name()
{
    chemkit::Fingerprint *fingerprint = chemkit::Fingerprint::create("pubchem");
    QVERIFY(fingerprint != 0);
    QCOMPARE(fingerprint->name(), std::string("pubchem"));
    delete fingerprint;
}

void PubChemTest::test_data()
{
    QTest::addColumn<QString>("cid");
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<QByteArray>("fingerprint_base64");

    QTest::newRow("acetylcholine") <<
        "187" <<
        "CC(=O)OCC[N+](C)(C)C" <<
        "C7H16NO2" <<
        QByteArray("AAADceBiMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHgAAAA"
                   "AAAADhgAYCCAMABAAIAACQCAAAAAAAAAAAAAAIAAACAAAAAAADAAAAAAAQ"
                   "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("adenosine") <<
        "60961" <<
        "C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)CO)O)O" <<
        "C10H13N5O4" <<
        QByteArray("AAADccBzuAAAAAAAAAAAAAAAAAAAAWJAAAAsAAAAAAAAAFgB+AAAHgAQCA"
                   "AACBzhlwYF8L9MFgCgAQZhZACAgC0REKABUCAoVBCDWAJAyEAeRAgPAALT"
                   "ACDwMAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("anthraquinone") <<
        "6780" <<
        "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2=O" <<
        "C14H8O2" <<
        QByteArray("AAADccBwMAAAAAAAAAAAAAAAAAAAAAAAAAAwYMAAAAAAAADBQAAAGgAAAA"
                   "AADASAmAAwAIAAAACIAqBSAAACAAAkAAAIiAEAAMgIIDKAFRCAIQAggAAI"
                   "iYcJiMCOgAAAAAAQAACAAAQAACAAAAAAAAAAAA==");
    QTest::newRow("caffeine") <<
        "2519" <<
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" <<
        "C8H10N4O2" <<
        QByteArray("AAADccBzsAAAAAAAAAAAAAAAAAAAAWAAAAAsAAAAAAAAAFgBgAAAHgAAAA"
                   "AACAgBlgQHsBcMEACoAQdxdACAgC0XEKABUAGoVECASAhASCAUAIgIByJA"
                   "AGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("cAMP") <<
        "6076" <<
        "C1C2C(C(C(O2)N3C=NC4=C3N=CN=C4N)O)OP(=O)(O1)O" <<
         "C10H12N5O6P" <<
         QByteArray("AAADccBzuAIAAAAAAAAAAAAAAAAAAWJAAAAsSAAAAAAAAFgB+AAAHgAQC"
                    "CAACBzhlwYF8L9MFxCgQQZhZICAgC0REKABUCAoVBCDWAJAyEAeRAgPAA"
                    "LTACDwMAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("diazepam") <<
        "3016" <<
        "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3" <<
        "C16H13ClN2O" <<
        QByteArray("AAADccB7IAAEAAAAAAAAAAAAAAAAAAAAAAAwYAAABAAAAAABQAAAHgIAAA"
                   "AADArBmCQywIMAAACoAyVyVACCAAAhBwAIiACoZpgIYKLBk5GUIAhghgDI"
                   "yAcQgMAOAABAIAACAAAAAIBAAAQAAAAAAAAAAA==");
    QTest::newRow("ethanol") <<
        "702" <<
        "CCO" <<
        "C2H6O" <<
        QByteArray("AAADcYBAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGgAACA"
                   "AAAACggAICAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAA"
                   "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("glucose") <<
        "5793" <<
        "C(C1C(C(C(C(O1)O)O)O)O)O" <<
        "C6H12O6" <<
        QByteArray("AAADccBgOAAAAAAAAAAAAAAAAAAAAAAAAAAkAAAAAAAAAAAAAAAAGgAACA"
                   "AACBSwgAMACAAABgAAAAAAAAAAAAAAAAAAAAAAAAAREAIAAAACQAAFAAAH"
                   "AAHAYAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("guanine") <<
        "764" <<
        "C1=NC2=C(N1)C(=O)N=C(N2)N" <<
        "C5H5N5O" <<
        QByteArray("AAADcYBjoAAAAAAAAAAAAAAAAAAAAWAAAAAgAAAAAAAAAEABgAAAHgAQAA"
                   "AACAgBlgQFsBbJkACoAQdxdACAgC2XEKABUYGoVECASAhASCAUAIAIAAJA"
                   "AGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("lysine") <<
        "866" <<
        "C(CCN)CC(C(=O)O)N" <<
        "C6H14N2O2" <<
        QByteArray("AAADccBjMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHgAQCA"
                   "AACCjBgAQACABAAgAIAACQCAAAAAAAAAAAAIGAAAACABIAgAAAQAAEEAAA"
                   "AAGYSAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");
    QTest::newRow("lsd") <<
        "5761" <<
        "CCN(CC)C(=O)C1CN(C2CC3=CNC4=CC=CC(=C34)C2=C1)C" <<
        "C20H25N3O" <<
        QByteArray("AAADceB7IAAAAAAAAAAAAAAAAAAAAWAAAAAwYIAAAAAAAFjB8AAAHgAQAA"
                   "AADSjBngQ+wPPJkACoAzV3VACCgCAxAiAI2aE4ZJgIIOrA0ZGEIAhglgDI"
                   "yAcQgMAOgAAAQAASAAAQAASAACQAAAAAAAAAAA==");
    QTest::newRow("octanitrocubane") <<
        "11762357" <<
         "C12(C3(C4(C1(C5(C2(C3(C45[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-])[N+]"
         "(=O)[O-])[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]" <<
         "C8N8O16" <<
         QByteArray("AAADcQBz/gAAAAAAAAAAAAAAAABgwAAAAAAAAAAAAAAAAAAAAAAADAAEA"
                    "AAADIgBAAAAAAAQQAABAAAAAwAAAAAAAAAgAAAAAAAAQAIAwAAAAAAAAA"
                    "AAAAEQgIAOgAAAAAAAAACQAQCACCQAQAAAAAAAAA==");
    QTest::newRow("thc") <<
        "16078" <<
        "CCCCCC1=CC2=C(C3C=C(CCC3C(O2)(C)C)C)C(=C1)O" <<
        "C21H30O2" <<
        QByteArray("AAADceB4MAAAAAAAAAAAAAAAAAAAAAAAAAA2QIAAAAAAAACRAAAAGgAACA"
                   "AADUSAmAAyBoAABgCAAiBCAAACCAAgIAAAiAAGCIgMJyKGMRqCeCClwBUI"
                   "uAeA4PwOwAABCAAIAACAAAIQABAAAAAAAAAAAA==");
    QTest::newRow("sertraline") <<
        "68617" <<
        "CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl" <<
        "C17H17Cl2N" <<
        QByteArray("AAADceB6AAAGAAAAAAAAAAAAAAAAAAAAAAAwYMAAAAAAAADBQAAAHAIQAA"
                   "AADSrBGCQyAILAAACAAiBCAACCAAAgBQAIisAIBogIICKBkxGEIAhgkAAI"
                   "iAcQgMAOhAAAIAAQAAQYAATAACQACAAAAAAAAA==");
    QTest::newRow("tnt") <<
        "8376" <<
        "CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]" <<
        "C7H5N3O6" <<
        QByteArray("AAADcYBjOAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAABAAAAHAAEAA"
                   "AADAiBGAAywIAQQACBAiRCQwCCAAAgAgAgiAAAZIoIICKA0dGAIABggAAI"
                   "yAcQgAAOCAAAAAQCAAAQAAAACAQAAAAAAAAAAA==");
}

void PubChemTest::test()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);
    QFETCH(QByteArray, fingerprint_base64);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());

    chemkit::Bitset fingerprint = molecule.fingerprint("pubchem");
    QCOMPARE(fingerprint.size(), size_t(881));

    QByteArray bytes = QByteArray::fromBase64(fingerprint_base64);

    // remove the first four bytes which contain the size
    bytes.remove(0, 4);

    boost::dynamic_bitset<unsigned char> expected;
    for(int i = 0; i < bytes.size(); i++){
        expected.append(reverse(static_cast<unsigned char>(bytes.at(i))));
    }

    compareFingerprints(fingerprint, expected);
}

QTEST_APPLESS_MAIN(PubChemTest)
