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

void StaticVectorTest::isZero()
{
    chemkit::StaticVector<float, 3> vector;
    QCOMPARE(vector.isZero(), true);

    vector.setValue(1, 1.0f);
    QCOMPARE(vector.isZero(), false);

    vector.setValue(1, 0.0f);
    QCOMPARE(vector.isZero(), true);
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
