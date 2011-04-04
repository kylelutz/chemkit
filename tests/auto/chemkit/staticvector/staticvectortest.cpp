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

#include "staticvectortest.h"

#include <chemkit/staticvector.h>

void StaticVectorTest::value()
{
    chemkit::StaticVector<float, 3> vector;
    vector.setValue(0, 3.0f);
    vector.setValue(1, 6.0f);
    vector.setValue(2, 9.0f);
    QCOMPARE(vector.value(0), 3.0f);
    QCOMPARE(vector.value(1), 6.0f);
    QCOMPARE(vector.value(2), 9.0f);
}

void StaticVectorTest::isNull()
{
    chemkit::StaticVector<float, 3> vector;
    QCOMPARE(vector.isNull(), true);

    vector.setValue(1, 1.0f);
    QCOMPARE(vector.isNull(), false);

    vector.setValue(1, 0.0f);
    QCOMPARE(vector.isNull(), true);
}

void StaticVectorTest::commaInitializer()
{
    chemkit::StaticVector<int, 4> vector;
    vector << 1, 2, 3, 4;
    QCOMPARE(vector[0], 1);
    QCOMPARE(vector[1], 2);
    QCOMPARE(vector[2], 3);
    QCOMPARE(vector[3], 4);

    vector << 4, 3, 2, 1;
    QCOMPARE(vector[0], 4);
    QCOMPARE(vector[1], 3);
    QCOMPARE(vector[2], 2);
    QCOMPARE(vector[3], 1);
}

void StaticVectorTest::dot()
{
    chemkit::StaticVector<int, 3> a;
    a << 1, 2, 3;

    chemkit::StaticVector<int, 3> b;
    b << 4, 5, 6;

    QCOMPARE(a.dot(b), 32);
    QCOMPARE(b.dot(a), 32);
    QCOMPARE(a * b, 32);
    QCOMPARE(b * a, 32);
}

void StaticVectorTest::cross()
{
    chemkit::StaticVector<int, 3> a;
    a << 1, 2, 3;

    chemkit::StaticVector<int, 3> b;
    b << 4, 5, 6;

    chemkit::StaticVector<int, 3> c;
    c << -3, 6, -3;

    QVERIFY(a.cross(b) == c);
    QVERIFY(b.cross(a) == -c);
    QVERIFY((a ^ b) == c);
    QVERIFY((b ^ a) == -c);
}

QTEST_APPLESS_MAIN(StaticVectorTest)
