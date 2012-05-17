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

#include "moriguchitest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void MoriguchiTest::initTestCase()
{
    // verify that the moriguchi-logp plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "moriguchi-logp") == 1);
}

void MoriguchiTest::logP_data()
{
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("logP");

    // data from examples in [Moriguchi 1992]
    QTest::newRow("halothane") << "C(C(F)(F)F)(Cl)Br" << "C2HBrClF3" << 2.60;
    QTest::newRow("ampicillin") << "O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H](c1ccccc1)N)[C@H]3SC2(C)C" << "C16H19N3O4S" << 1.06;
    QTest::newRow("valerolactone") << "CC1CCC(=O)O1" << "C5H8O2" << 0.60;
    QTest::newRow("oxazepam") << "C1=CC=C(C=C1)C2=NC(C(=O)NC3=C2C=C(C=C3)Cl)O" << "C15H11ClN2O2" << 3.12;

    // data from Table 2 in [Moriguchi 1994]
    QTest::newRow("atropine") << "CN3[C@H]1CC[C@@H]3C[C@@H](C1)OC(=O)C(CO)c2ccccc2" << "C17H23NO3" << 2.21;
    QTest::newRow("chloramphenicol") << "c1cc(ccc1[C@H]([C@@H](CO)NC(=O)C(Cl)Cl)O)[N+](=O)[O-]" << "C11H12Cl2N2O5" << 1.23;
    //QTest::newRow("chlorothiazide") << "O=S(=O)(c1c(Cl)cc2c(c1)S(=O)(=O)/N=C\\N2)N" << "C7H6ClN3O4S2" << -0.36;
    //QTest::newRow("chlorpromazine") << "CN(C)CCCN1c2ccccc2Sc3c1cc(cc3)Cl" << "C17H19ClN2S" << 3.77;
    QTest::newRow("cimetidine") << "N#CN\\C(=N/C)NCCSCc1ncnc1C" << "C10H16N6S" << 0.82;
    QTest::newRow("diazepam") << "CN1c2ccc(cc2C(=NCC1=O)c3ccccc3)Cl" << "C16H13ClN2O" << 3.36;
    QTest::newRow("diltiazem") << "O=C2N(c3c(S[C@@H](c1ccc(OC)cc1)[C@H]2OC(=O)C)cccc3)CCN(C)C" << "C22H26N2O4S" << 2.67;
    QTest::newRow("diphenhydramine") << "O(CCN(C)C)C(c1ccccc1)c2ccccc2" << "C17H21NO" << 3.26;
    QTest::newRow("disopyramide") << "O=C(N)C(c1ncccc1)(c2ccccc2)CCN(C(C)C)C(C)C" << "C21H29N3O" << 2.47;
    QTest::newRow("flufenamic acid") << "FC(F)(F)c1cc(ccc1)Nc2ccccc2C(=O)O" << "C14H10F3NO2" << 3.86;
    QTest::newRow("haloperidol") << "c1cc(ccc1C(=O)CCCN2CCC(CC2)(c3ccc(cc3)Cl)O)F" << "C21H23ClFNO2" << 4.01;
    QTest::newRow("imipramine") << "c1cc3c(cc1)CCc2c(cccc2)N3CCCN(C)C" << "C19H24N2" << 3.88;
    QTest::newRow("lidocaine") << "O=C(Nc1c(cccc1C)C)CN(CC)CC" << "C14H22N2O" << 2.52;
    QTest::newRow("phenobarbital") << "O=C1NC(=O)NC(=O)C1(c2ccccc2)CC" << "C12H12N2O3" << 0.78;
    QTest::newRow("phenytoin") << "O=C2NC(=O)NC2(c1ccccc1)c3ccccc3" << "C15H12N2O2" << 1.80;
    QTest::newRow("procainamide") << "O=C(c1ccc(N)cc1)NCCN(CC)CC" << "C13H21N3O" << 1.72;
    QTest::newRow("propafenone") << "O=C(c1ccccc1OCC(O)CNCCC)CCc2ccccc2" << "C21H27NO3" << 3.00;
    QTest::newRow("propranolol") << "CC(C)NCC(COc1cccc2c1cccc2)O" << "C16H21NO2" << 2.53;
    QTest::newRow("tetracaine") << "O=C(OCCN(C)C)c1ccc(NCCCC)cc1" << "C15H24N2O2" << 2.64;
    //QTest::newRow("trimethoprim") << "COc1cc(cc(c1OC)OC)Cc2cnc(nc2N)N" << "C14H18N4O3" << 1.26;
    QTest::newRow("verapamil") << "N#CC(c1cc(OC)c(OC)cc1)(CCCN(CCc2ccc(OC)c(OC)c2)C)C(C)C" << "C27H38N2O4" << 3.23;
}

void MoriguchiTest::logP()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);
    QFETCH(double, logP);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());

    double tolerance = 0.1;
    double actual = molecule.descriptor("moriguchi-logp").toDouble();
    bool equal = std::abs(logP - actual) < tolerance;
    if(!equal){
        qDebug() << "Actual: " << actual;
        qDebug() << "Expected: " << logP;
    }
    QVERIFY(equal);
}

QTEST_APPLESS_MAIN(MoriguchiTest)
