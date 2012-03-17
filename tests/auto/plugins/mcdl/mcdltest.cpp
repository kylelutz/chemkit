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

#include "mcdltest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void McdlTest::initTestCase()
{
    // verify that the mcdl plugin registered itself correctly
    QVERIFY(boost::count(chemkit::LineFormat::formats(), "mcdl") == 1);
}

void McdlTest::read_data()
{
    QTest::addColumn<QString>("mcdlString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<int>("atomCount");
    QTest::addColumn<int>("bondCount");
    QTest::addColumn<int>("ringCount");

    QTest::newRow("acetylChloride") << "CHHH;COCl[2]" << "C2H3ClO" << 7 << 6 << 0;
    QTest::newRow("adenine") << "3C;2CH;3N;NH;NHH[2,3,6;7,9;8,10;7,8;6,9]" << "C5H5N5" << 15 << 16 << 2;
    QTest::newRow("asprin") << "2C;4CH;CHHH;2CO;O;OH[2,3,8;4,10;5;6;6;;9;11;10]" << "C9H8O4" << 21 << 21 << 1;
    QTest::newRow("bromouracil") << "CBr;CH;2CO;2NH[2,3;5;6;5,6]" << "C4H3BrN2O2" << 12 << 12 << 1;
    QTest::newRow("caffeine") << "2C;CH;3CHHH;2CO;4N[2,7,9;10,11;9,10;9;11;12;12;11,12]" << "C8H10N4O2" << 24 << 25 << 2;
    QTest::newRow("ethanol") << "CHH;CHHH;OH[2,3]" << "C2H6O" << 9 << 8 << 0;
    QTest::newRow("guanine") << "3C;CH;CO;2N;2NH;NHH[2,5,7;6,8;6,9,10;7,8;9]" << "C5H5N5O" << 16 << 17 << 2;
    QTest::newRow("indole") << "C;7CH;N[2,3,9;4;5;6;7;8;9;9]" << "C8H7N" << 16 << 17 << 2;
    QTest::newRow("uridine") << "6CH;CHH;2CO;N;NH;O;3OH[2,3,13;5,14;7,12;6,8;10,12;10;15;11;10,11]" << "C9H12N2O6" << 29 << 30 << 2;
    QTest::newRow("water") << "OHH[]" << "H2O" << 3 << 2 << 0;
}

void McdlTest::read()
{
    chemkit::LineFormat *mcdlFormat = chemkit::LineFormat::create("mcdl");
    if(!mcdlFormat){
        return;
    }

    QFETCH(QString, mcdlString);
    QFETCH(QString, formulaString);
    QFETCH(int, atomCount);
    QFETCH(int, bondCount);
    QFETCH(int, ringCount);

    QByteArray mcdl = mcdlString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule *molecule = mcdlFormat->read(mcdl.constData());
    if(!molecule)
        qDebug() << mcdlFormat->errorString().c_str();
    QVERIFY(molecule != 0);

    QCOMPARE(molecule->formula().c_str(), formula.constData());
    QCOMPARE(molecule->atomCount(), size_t(atomCount));
    QCOMPARE(molecule->bondCount(), size_t(bondCount));
    QCOMPARE(molecule->ringCount(), size_t(ringCount));

    delete molecule;
    delete mcdlFormat;
}

QTEST_APPLESS_MAIN(McdlTest)
