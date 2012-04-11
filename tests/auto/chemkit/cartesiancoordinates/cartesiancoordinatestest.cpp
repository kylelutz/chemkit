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

#include "cartesiancoordinatestest.h"

#include <chemkit/atom.h>
#include <chemkit/vector3.h>
#include <chemkit/molecule.h>
#include <chemkit/cartesiancoordinates.h>

void CartesianCoordinatesTest::setPosition()
{
    chemkit::CartesianCoordinates matrix(5);
    QCOMPARE(matrix.position(0), chemkit::Point3(0, 0, 0));

    matrix.setPosition(1, chemkit::Point3(1, 2, 3));
    QCOMPARE(matrix.position(1), chemkit::Point3(1, 2, 3));

    matrix.setPosition(2, -5, 8, 0.5);
    QCOMPARE(matrix.position(2), chemkit::Point3(-5, 8, 0.5));
}

void CartesianCoordinatesTest::append()
{
    chemkit::CartesianCoordinates matrix;
    QCOMPARE(matrix.size(), size_t(0));

    matrix.append(1, 2, 3);
    QCOMPARE(matrix.size(), size_t(1));
    QCOMPARE(matrix.position(0), chemkit::Point3(1, 2, 3));

    matrix.append(4, 5, 6);
    QCOMPARE(matrix.size(), size_t(2));
    QCOMPARE(matrix.position(0), chemkit::Point3(1, 2, 3));
    QCOMPARE(matrix.position(1), chemkit::Point3(4, 5, 6));
}

void CartesianCoordinatesTest::insert()
{
    chemkit::CartesianCoordinates matrix(3);
    matrix.setPosition(0, 1.0, 2.0, 3.0);
    matrix.setPosition(1, 4.0, 5.0, 6.0);
    matrix.setPosition(2, 7.0, 8.0, 9.0);

    matrix.insert(3, -1.0, -2.0, -3.0);
    QCOMPARE(matrix.size(), size_t(4));
    QCOMPARE(matrix.position(3), chemkit::Point3(-1.0, -2.0, -3.0));
    QCOMPARE(matrix.position(2), chemkit::Point3(7.0, 8.0, 9.0));

    matrix.insert(1, 0.5, 1.5, 2.5);
    QCOMPARE(matrix.size(), size_t(5));
    QCOMPARE(matrix.position(0), chemkit::Point3(1.0, 2.0, 3.0));
    QCOMPARE(matrix.position(1), chemkit::Point3(0.5, 1.5, 2.5));
    QCOMPARE(matrix.position(2), chemkit::Point3(4.0, 5.0, 6.0));
    QCOMPARE(matrix.position(3), chemkit::Point3(7.0, 8.0, 9.0));
}

void CartesianCoordinatesTest::remove()
{
    chemkit::CartesianCoordinates matrix(3);
    matrix.setPosition(0, 1.0, 2.0, 3.0);
    matrix.setPosition(1, 4.0, 5.0, 6.0);
    matrix.setPosition(2, 7.0, 8.0, 9.0);

    matrix.remove(0);
    QCOMPARE(matrix.size(), size_t(2));
    QCOMPARE(matrix.position(0), chemkit::Point3(4.0, 5.0, 6.0));
    QCOMPARE(matrix.position(1), chemkit::Point3(7.0, 8.0, 9.0));

    matrix.remove(1);
    QCOMPARE(matrix.size(), size_t(1));
    QCOMPARE(matrix.position(0), chemkit::Point3(4.0, 5.0, 6.0));

    matrix.remove(0);
    QCOMPARE(matrix.size(), size_t(0));
}

void CartesianCoordinatesTest::center()
{
    chemkit::CartesianCoordinates matrix;
    QCOMPARE(matrix.center(), chemkit::Point3(0, 0, 0));

    matrix = chemkit::CartesianCoordinates(2);
    matrix.setPosition(0, chemkit::Point3(0, 0, 0));
    matrix.setPosition(1, chemkit::Point3(0, 5, 0));
    QCOMPARE(matrix.center(), chemkit::Point3(0, 2.5, 0));
}

void CartesianCoordinatesTest::multiply()
{
    chemkit::CartesianCoordinates a(7);
    a.setPosition(0, chemkit::Point3(5, 5, 5));
    a.setPosition(1, chemkit::Point3(8, -2, 1.5));
    a.setPosition(2, chemkit::Point3(0, 4, 1));
    a.setPosition(3, chemkit::Point3(-1, -3, 8));
    a.setPosition(4, chemkit::Point3(2, 10, 15));
    a.setPosition(5, chemkit::Point3(-1, 2.5, 3));
    a.setPosition(6, chemkit::Point3(0, -9, 11.75));

    chemkit::CartesianCoordinates b(7);
    b.setPosition(0, chemkit::Point3(19, 12, 1));
    b.setPosition(1, chemkit::Point3(0, 0, 0));
    b.setPosition(2, chemkit::Point3(-8, -9, 13));
    b.setPosition(3, chemkit::Point3(9, 8, 1.4));
    b.setPosition(4, chemkit::Point3(6.7, -3, -4.2));
    b.setPosition(5, chemkit::Point3(0, 8, 9));
    b.setPosition(6, chemkit::Point3(-2.5, 1.5, 0));

    Eigen::Matrix<chemkit::Real, 3, 3> product = a.multiply(&b);
    QCOMPARE(product(0, 0), chemkit::Real(99.4));
    QCOMPARE(product(0, 1), chemkit::Real(38.0));
    QCOMPARE(product(0, 2), chemkit::Real(-13.8));
    QCOMPARE(product(1, 0), chemkit::Real(125.5));
    QCOMPARE(product(1, 1), chemkit::Real(-23.5));
    QCOMPARE(product(1, 2), chemkit::Real(33.3));
    QCOMPARE(product(2, 0), chemkit::Real(230.125));
    QCOMPARE(product(2, 1), chemkit::Real(111.625));
    QCOMPARE(product(2, 2), chemkit::Real(-6.8));
}

void CartesianCoordinatesTest::distanceMatrix()
{
    chemkit::CartesianCoordinates coordinates(4);
    coordinates.setPosition(0, chemkit::Point3(1, 0, 0));
    coordinates.setPosition(1, chemkit::Point3(2, 0, 0));
    coordinates.setPosition(2, chemkit::Point3(0, 5, 0));
    coordinates.setPosition(3, chemkit::Point3(10, 5, 2));

    chemkit::Matrix distances = coordinates.distanceMatrix();
    QVERIFY(distances.rows() == 4);
    QVERIFY(distances.cols() == 4);
    QCOMPARE(qRound(distances(0, 0)), 0);
    QCOMPARE(qRound(distances(0, 1)), 1);
    QCOMPARE(qRound(distances(1, 0)), 1);
}

void CartesianCoordinatesTest::rotate()
{
    chemkit::CartesianCoordinates coordinates(3);
    coordinates.setPosition(0, chemkit::Point3(0, -1, 0));
    coordinates.setPosition(1, chemkit::Point3(0, 0, 0));
    coordinates.setPosition(2, chemkit::Point3(1, 0, 0));

    coordinates.rotate(chemkit::Vector3::UnitZ(), 90);
    QVERIFY(coordinates.position(0).isApprox(chemkit::Point3(1, 0, 0)));
    QVERIFY(coordinates.position(1).isApprox(chemkit::Point3(0, 0, 0)));
    QVERIFY(coordinates.position(2).isApprox(chemkit::Point3(0, 1, 0)));
}

QTEST_APPLESS_MAIN(CartesianCoordinatesTest)
