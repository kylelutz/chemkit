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

#include "point3test.h"

#include <chemkit/point3.h>

void Point3Test::basic()
{
    chemkit::Point3 point;
    QCOMPARE(point.x(), chemkit::Float(0.0));
    QCOMPARE(point.y(), chemkit::Float(0.0));
    QCOMPARE(point.z(), chemkit::Float(0.0));

    point = chemkit::Point3(1, 2, 3);
    QCOMPARE(point.x(), chemkit::Float(1.0));
    QCOMPARE(point.y(), chemkit::Float(2.0));
    QCOMPARE(point.z(), chemkit::Float(3.0));
}

void Point3Test::moveBy()
{
}

void Point3Test::movedBy()
{
}

void Point3Test::isZero()
{
    chemkit::Point3 point(0, 0, 0);
    QCOMPARE(point.isZero(), true);

    point = chemkit::Point3(1, 2, 3);
    QCOMPARE(point.isZero(), false);

    point = chemkit::Point3(0, 0, 0);
    QCOMPARE(point.isZero(), true);
}

void Point3Test::distance()
{
    chemkit::Point3 a, b;
    QCOMPARE(a.distance(b), chemkit::Float(0.0));

    a = chemkit::Point3(1, 0, 0);
    b = chemkit::Point3(3, 0, 0);
    QCOMPARE(a.distance(b), chemkit::Float(2.0));
}

void Point3Test::midpoint()
{
    chemkit::Point3 a(1, 0, 0);
    chemkit::Point3 b(3, 0, 0);

    chemkit::Point3 midpoint = a.midpoint(b);
    QCOMPARE(midpoint, chemkit::Point3(2, 0, 0));
}

QTEST_APPLESS_MAIN(Point3Test)
