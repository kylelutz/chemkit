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
