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

#include "mcdltest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void McdlTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::LineFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "mcdl") != formats.end());
}

void McdlTest::read_data()
{
    QTest::addColumn<QString>("mcdl");
    QTest::addColumn<QString>("formula");
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

    QFETCH(QString, mcdl);
    QFETCH(QString, formula);
    QFETCH(int, atomCount);
    QFETCH(int, bondCount);
    QFETCH(int, ringCount);

    chemkit::Molecule molecule;
    bool ok = mcdlFormat->read(mcdl.toStdString(), &molecule);
    if(!ok)
        qDebug() << mcdlFormat->errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(molecule.formula(), formula.toStdString());
    QCOMPARE(molecule.atomCount(), atomCount);
    QCOMPARE(molecule.bondCount(), bondCount);
    QCOMPARE(molecule.ringCount(), ringCount);

    delete mcdlFormat;
}

QTEST_APPLESS_MAIN(McdlTest)
