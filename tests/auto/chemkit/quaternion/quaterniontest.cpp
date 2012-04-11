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
    chemkit::Quaternion q(4, 1, 2, 3);
    QCOMPARE(q.x(), chemkit::Real(1.0));
    QCOMPARE(q.y(), chemkit::Real(2.0));
    QCOMPARE(q.z(), chemkit::Real(3.0));
    QCOMPARE(q.w(), chemkit::Real(4.0));

    q = chemkit::Quaternion(8, 5, 6, 7);
    QCOMPARE(q.x(), chemkit::Real(5.0));
    QCOMPARE(q.y(), chemkit::Real(6.0));
    QCOMPARE(q.z(), chemkit::Real(7.0));
    QCOMPARE(q.w(), chemkit::Real(8.0));
}

void QuaternionTest::multiply()
{
    chemkit::Quaternion a(5, 2, 3, 4);
    chemkit::Quaternion b(7, 4, 5, 6);
    chemkit::Quaternion c = a * b;
    QCOMPARE(qRound(c.x()), 32);
    QCOMPARE(qRound(c.y()), 50);
    QCOMPARE(qRound(c.z()), 56);
    QCOMPARE(qRound(c.w()), -12);

    c = b * a;
    QCOMPARE(qRound(c.x()), 36);
    QCOMPARE(qRound(c.y()), 42);
    QCOMPARE(qRound(c.z()), 60);
    QCOMPARE(qRound(c.w()), -12);
}

void QuaternionTest::conjugate()
{
    chemkit::Quaternion q(4, 1, 2, 3);
    chemkit::Quaternion p = q.conjugate();
    QCOMPARE(qRound(p.x()), -1);
    QCOMPARE(qRound(p.y()), -2);
    QCOMPARE(qRound(p.z()), -3);
    QCOMPARE(qRound(p.w()), 4);

    q = chemkit::Quaternion(-8, -2, -4, -6);
    p = q.conjugate();
    QCOMPARE(qRound(p.x()), 2);
    QCOMPARE(qRound(p.y()), 4);
    QCOMPARE(qRound(p.z()), 6);
    QCOMPARE(qRound(p.w()), -8);
}

QTEST_APPLESS_MAIN(QuaternionTest)
