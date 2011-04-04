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

#include "graphicstransformtest.h"

#include <chemkit/graphicstransform.h>

void GraphicsTransformTest::data()
{
}

void GraphicsTransformTest::multiplyPoint()
{
    chemkit::Point3f point(1, 2, 3);
    chemkit::GraphicsTransform identity = chemkit::GraphicsTransform::identity();

    point = identity.multiply(point);
    QCOMPARE(point, chemkit::Point3f(1, 2, 3));

    chemkit::GraphicsTransform rotation = chemkit::GraphicsTransform::rotation(chemkit::Vector3f::X(), 180);
    point = rotation.multiply(point);
    QCOMPARE(point, chemkit::Point3f(1, -2, -3));
}

void GraphicsTransformTest::multiplyVector()
{
    chemkit::Vector3f vector(1, 2, 3);
    chemkit::GraphicsTransform identity = chemkit::GraphicsTransform::identity();

    vector = identity.multiply(vector);
    QCOMPARE(vector, chemkit::Vector3f(1, 2, 3));

    // translation matrix should have no effect on vectors
    chemkit::Vector3f translationVector(5, 5, 5);
    chemkit::GraphicsTransform translation = chemkit::GraphicsTransform::translation(translationVector);
    vector = translation.multiply(vector);
    QCOMPARE(vector, chemkit::Vector3f(1, 2, 3));
}

void GraphicsTransformTest::multiplyTransform()
{
    // identity transform
    chemkit::GraphicsTransform I = chemkit::GraphicsTransform::identity();

    // fill transform A with the values from 1 to 16
    float A_data[] = {1, 2, 3, 4,
                      5, 6, 7, 8,
                      9, 10, 11, 12,
                      13, 14, 15, 16};
    chemkit::GraphicsTransform A = chemkit::GraphicsTransform(A_data);

    // multiply by identity matrix
    A = I * A;

    // ensure values remain the same
    QCOMPARE(qRound(A(0, 0)), 1);
    QCOMPARE(qRound(A(1, 2)), 7);
    QCOMPARE(qRound(A(3, 3)), 16);

    // multiply A by itself
    A = A * A;
    QCOMPARE(qRound(A(0, 0)), 90);
    QCOMPARE(qRound(A(1, 2)), 254);
    QCOMPARE(qRound(A(3, 3)), 600);
}

void GraphicsTransformTest::inverseMultiplyPoint()
{
}

void GraphicsTransformTest::inverseMultiplyVector()
{
}

void GraphicsTransformTest::translation()
{
    chemkit::Vector3f translationVector(5, 4, 3);
    chemkit::GraphicsTransform transform = chemkit::GraphicsTransform::translation(translationVector);
    QCOMPARE(qRound(transform(0, 3)), 5);
    QCOMPARE(qRound(transform(1, 3)), 4);
    QCOMPARE(qRound(transform(2, 3)), 3);
    QCOMPARE(qRound(transform(3, 3)), 1);

    chemkit::Point3f point(0, 0, 0);
    QCOMPARE(point, chemkit::Point3f(0, 0, 0));

    point = transform * point;
    QCOMPARE(point, chemkit::Point3f(5, 4, 3));

    point = transform * point;
    QCOMPARE(point, chemkit::Point3f(10, 8, 6));

    point = transform * point;
    QCOMPARE(point, chemkit::Point3f(15, 12, 9));
}

QTEST_APPLESS_MAIN(GraphicsTransformTest)
