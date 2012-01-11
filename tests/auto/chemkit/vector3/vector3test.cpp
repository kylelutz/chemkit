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

#include "vector3test.h"

#include <chemkit/vector3.h>

void Vector3Test::constructor()
{
    chemkit::Vector3 empty(0, 0, 0);
    QCOMPARE(qRound(empty.x()), 0);
    QCOMPARE(qRound(empty.y()), 0);
    QCOMPARE(qRound(empty.z()), 0);

    chemkit::Vector3 vector123(1, 2, 3);
    QCOMPARE(qRound(vector123.x()), 1);
    QCOMPARE(qRound(vector123.y()), 2);
    QCOMPARE(qRound(vector123.z()), 3);
}

void Vector3Test::value()
{
    chemkit::Vector3 vector(0, 0, 0);
    QCOMPARE(qRound(vector[0]), 0);
    QCOMPARE(qRound(vector[1]), 0);
    QCOMPARE(qRound(vector[2]), 0);

    vector[0] = 1;
    QCOMPARE(qRound(vector[0]), 1);

    vector[1] = 3;
    QCOMPARE(qRound(vector[1]), 3);

    vector[2] = -5;
    QCOMPARE(qRound(vector[2]), -5);

    vector[0] = vector[1];
    QCOMPARE(qRound(vector[0]), 3);
}

void Vector3Test::norm()
{
    chemkit::Vector3 vector(0, 0, 0);
    QCOMPARE(qRound(vector.norm()), 0);

    vector = chemkit::Vector3(1, 0, 0);
    QCOMPARE(qRound(vector.x()), 1);
    QCOMPARE(qRound(vector.y()), 0);
    QCOMPARE(qRound(vector.z()), 0);
    QCOMPARE(qRound(vector.norm()), 1);
}

void Vector3Test::isZero()
{
    chemkit::Vector3 vector(0, 0, 0);
    QCOMPARE(vector.isZero(), true);

    vector = chemkit::Vector3(1, 2, 3);
    QCOMPARE(vector.isZero(), false);
}

void Vector3Test::normalize()
{
    chemkit::Vector3 vector(2, 0, 0);
    QCOMPARE(qRound(vector.norm()), 2);

    vector.normalize();
    QCOMPARE(qRound(vector.norm()), 1);
    QCOMPARE(qRound(vector.x()), 1);

    chemkit::Vector3 nullVector(0, 0, 0);
    QCOMPARE(qRound(nullVector.norm()), 0);
    nullVector.normalize();
    QCOMPARE(qRound(nullVector.norm()), 0);
}

QTEST_APPLESS_MAIN(Vector3Test)
