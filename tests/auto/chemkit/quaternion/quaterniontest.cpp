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

#include "quaterniontest.h"

#include <chemkit/quaternion.h>

void QuaternionTest::basic()
{
    chemkit::Quaternion q(1, 2, 3, 4);
    QCOMPARE(q.x(), 1.0);
    QCOMPARE(q.y(), 2.0);
    QCOMPARE(q.z(), 3.0);
    QCOMPARE(q.r(), 4.0);

    q = chemkit::Quaternion(5, 6, 7, 8);
    QCOMPARE(q.x(), 5.0);
    QCOMPARE(q.y(), 6.0);
    QCOMPARE(q.z(), 7.0);
    QCOMPARE(q.r(), 8.0);
}

void QuaternionTest::add()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.add(b) == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(a + b == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(b.add(a) == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(b + a == chemkit::Quaternion(6, 8, 10, 12));
}

void QuaternionTest::subtract()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.subtract(b) == chemkit::Quaternion(-2, -2, -2, -2));
    QVERIFY(a - b == chemkit::Quaternion(-2, -2, -2, -2));
    QVERIFY(b.subtract(a) == chemkit::Quaternion(2, 2, 2, 2));
    QVERIFY(b - a == chemkit::Quaternion(2, 2, 2, 2));
}

void QuaternionTest::multiply()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.multiply(b) == chemkit::Quaternion(32, 50, 56, -12));
    QVERIFY(a * b == chemkit::Quaternion(32, 50, 56, -12));
    QVERIFY(b.multiply(a) == chemkit::Quaternion(36, 42, 60, -12));
    QVERIFY(b * a == chemkit::Quaternion(36, 42, 60, -12));
}

void QuaternionTest::conjugate()
{
    chemkit::Quaternion q(1, 2, 3, 4);
    QVERIFY(q.conjugate() == chemkit::Quaternion(-1, -2, -3, 4));

    q = chemkit::Quaternion(-2, -4, -6, -8);
    QVERIFY(q.conjugate() == chemkit::Quaternion(2, 4, 6, -8));
}

QTEST_APPLESS_MAIN(QuaternionTest)
