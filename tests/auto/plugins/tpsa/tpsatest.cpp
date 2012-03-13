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

#include "tpsatest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void TpsaTest::initTestCase()
{
    // verify that the tpsa plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "tpsa") == 1);
}

void TpsaTest::test_data()
{
    QTest::addColumn<QString>("smiles");
    QTest::addColumn<double>("tpsa");

    // data from table 3 in [Ertl 2000]
    QTest::newRow("metoprolol") << "O(c1ccc(cc1)CCOC)CC(O)CNC(C)C" << 50.7;
    QTest::newRow("nordiazepam") << "c1ccc(cc1)C2=NCC(=O)Nc3c2cc(cc3)Cl" << 41.5;
    QTest::newRow("diazepam") << "CN1c2ccc(cc2C(=NCC1=O)c3ccccc3)Cl" << 32.7;
    QTest::newRow("oxprenolol") << "CC(C)NCC(COC1=CC=CC=C1OCC=C)O" << 50.7;
    QTest::newRow("phenazone") << "CC1=CC(=O)N(N1C)C2=CC=CC=C2" << 26.9;
    QTest::newRow("oxazepam") << "C1=CC=C(C=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O" << 61.7;
    QTest::newRow("alprenolol") << "CC(C)NCC(COC1=CC=CC=C1CC=C)O" << 41.9;
    QTest::newRow("practolol") << "CC(C)NCC(COC1=CC=C(C=C1)NC(=O)C)O" << 70.6;
    QTest::newRow("pindolol") << "CC(C)NCC(COC1=CC=CC2=C1C=CN2)O" << 57.3;
    QTest::newRow("ciprofloxacin") << "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O" << 72.9; // different due to aromaticity
    QTest::newRow("metolazone") << "CC1NC2=CC(=C(C=C2C(=O)N1C3=CC=CC=C3C)S(=O)(=O)N)Cl" << 92.5;
    QTest::newRow("tranexamic acid") << "C1CC(CCC1CN)C(=O)O" << 63.3;
    QTest::newRow("atenolol") << "CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O" << 84.6;
    QTest::newRow("sulpiride") << "CCN1CCCC1CNC(=O)C2=C(C=CC(=C2)S(=O)(=O)N)OC" << 101.7;
    QTest::newRow("mannitol") << "C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O" << 121.4;
    QTest::newRow("foscarnet") << "C(=O)(O)P(=O)(O)O" << 104.64; // table 3 value (94.8) does not include P
    QTest::newRow("sulfasalazine") << "c1ccnc(NS(=O)(c2ccc(N=Nc3ccc(c(C(O)=O)c3)O)cc2)=O)c1" << 141.3;
    QTest::newRow("olsalazine") << "c1cc(O)c(C(O)=O)cc1N=Nc2ccc(c(C(O)=O)c2)O" << 139.8;
    QTest::newRow("lactulose") << "C(O)C1C(O)C(O)C(O)C(OC(C(O)C(O)CO)C(CO)=O)O1" << 197.4;
    QTest::newRow("raffinose") << "C(C1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3(C(C(C(O3)CO)O)O)CO)O)O)O)O)O)O)O" << 268.7;

    // data from pubchem
    QTest::newRow("ethanol") << "CCO" << 20.2;
    QTest::newRow("formaldehyde") << "C=O" << 17.1;
    QTest::newRow("glycine") << "C(C(=O)O)N" << 63.3;
    QTest::newRow("alanine") << "CC(C(=O)O)N" << 63.3;
    QTest::newRow("asparagine") << "C(C(C(=O)O)N)C(=O)N" << 106.0;
    QTest::newRow("glutamic acid") << "C(CC(=O)O)C(C(=O)O)N" << 101.0;
    QTest::newRow("paracetamol") << "CC(=O)NC1=CC=C(C=C1)O" << 49.3;
    QTest::newRow("uracil") << "C1=CNC(=O)NC1=O" << 65.72; // different due to aromaticity
    QTest::newRow("adenosine") << "C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)CO)O)O" << 139.54;
    QTest::newRow("cAMP") << "C1C2C(C(C(O2)N3C=NC4=C3N=CN=C4N)O)OP(=O)(O1)O" << 164.65; // different due to not including P
    QTest::newRow("trinitrotoluene") << "CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]" << 137.46;
}

void TpsaTest::test()
{
    QFETCH(QString, smiles);
    QFETCH(double, tpsa);

    chemkit::Molecule molecule(smiles.toStdString(), "smiles");
    QVERIFY(!molecule.isEmpty());

    const double tolerance = 0.5;
    double actual = molecule.descriptor("tpsa").toDouble();

    bool equal = std::abs(actual - tpsa) < tolerance;
    if(!equal){
        qDebug() << "Actual: " << actual;
        qDebug() << "Expected: " << tpsa;
        QVERIFY(equal);
    }
}

QTEST_APPLESS_MAIN(TpsaTest)
