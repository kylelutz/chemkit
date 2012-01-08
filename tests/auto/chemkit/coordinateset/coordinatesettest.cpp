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

#include "coordinatesettest.h"

#include <chemkit/coordinateset.h>
#include <chemkit/diagramcoordinates.h>
#include <chemkit/internalcoordinates.h>
#include <chemkit/cartesiancoordinates.h>

void CoordinateSetTest::type()
{
    chemkit::CoordinateSet coordinates;
    QVERIFY(coordinates.type() == chemkit::CoordinateSet::None);

    chemkit::CoordinateSet cartesianCoordinates(new chemkit::CartesianCoordinates);
    QVERIFY(cartesianCoordinates.type() == chemkit::CoordinateSet::Cartesian);

    chemkit::CoordinateSet internalCoordinates(new chemkit::InternalCoordinates);
    QVERIFY(internalCoordinates.type() == chemkit::CoordinateSet::Internal);

    chemkit::CoordinateSet diagramCoordinates(new chemkit::DiagramCoordinates);
    QVERIFY(diagramCoordinates.type() == chemkit::CoordinateSet::Diagram);
}

void CoordinateSetTest::size()
{
    chemkit::CoordinateSet coordinates;
    QCOMPARE(coordinates.size(), size_t(0));

    coordinates = chemkit::CoordinateSet(new chemkit::CartesianCoordinates(10));
    QCOMPARE(coordinates.size(), size_t(10));
}

void CoordinateSetTest::isEmpty()
{
    chemkit::CoordinateSet coordinates;
    QCOMPARE(coordinates.isEmpty(), true);

    coordinates = chemkit::CoordinateSet(new chemkit::CartesianCoordinates(10));
    QCOMPARE(coordinates.isEmpty(), false);
}

void CoordinateSetTest::cartesianCoordinates()
{
    chemkit::CoordinateSet coordinates;
    QVERIFY(coordinates.cartesianCoordinates() == 0);
}

void CoordinateSetTest::internalCoordinates()
{
    chemkit::CoordinateSet coordinates;
    QVERIFY(coordinates.internalCoordinates() == 0);
}

void CoordinateSetTest::diagramCoordinates()
{
    chemkit::CoordinateSet coordinates;
    QVERIFY(coordinates.diagramCoordinates() == 0);
}

void CoordinateSetTest::position()
{
    chemkit::CoordinateSet coordinates;
    coordinates.setCoordinates(new chemkit::CartesianCoordinates(3));
    coordinates.cartesianCoordinates()->setPosition(0, chemkit::Point3(1, 2, 3));
    coordinates.cartesianCoordinates()->setPosition(1, chemkit::Point3(4, 5, 6));
    coordinates.cartesianCoordinates()->setPosition(2, chemkit::Point3(7, 8, 9));
    QCOMPARE(coordinates.position(0), chemkit::Point3(1, 2, 3));
    QCOMPARE(coordinates.position(1), chemkit::Point3(4, 5, 6));
    QCOMPARE(coordinates.position(2), chemkit::Point3(7, 8, 9));

    coordinates.setCoordinates(new chemkit::DiagramCoordinates(3));
    coordinates.diagramCoordinates()->setPosition(0, chemkit::Point2f(10, 15));
    coordinates.diagramCoordinates()->setPosition(1, chemkit::Point2f(30, 45));
    coordinates.diagramCoordinates()->setPosition(2, chemkit::Point2f(60, 75));
    QCOMPARE(coordinates.position(0), chemkit::Point3(10, 15, 0));
    QCOMPARE(coordinates.position(1), chemkit::Point3(30, 45, 0));
    QCOMPARE(coordinates.position(2), chemkit::Point3(60, 75, 0));
}

QTEST_APPLESS_MAIN(CoordinateSetTest)
