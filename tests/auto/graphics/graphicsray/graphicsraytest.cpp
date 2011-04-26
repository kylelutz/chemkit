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

#include "graphicsraytest.h"

#include <chemkit/graphicsray.h>

void GraphicsRayTest::basic()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.origin(), chemkit::Point3f(0, 0, 0));
    QCOMPARE(ray.direction(), chemkit::Vector3f(0, 0, -1));

    ray = chemkit::GraphicsRay(chemkit::Point3f(0, 0, 0), chemkit::Vector3f(0, 1, 0));
    QCOMPARE(ray.origin(), chemkit::Point3f(0, 0, 0));
    QCOMPARE(ray.direction(), chemkit::Vector3f(0, 1, 0));
}

void GraphicsRayTest::setOrigin()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.origin(), chemkit::Point3f(0, 0, 0));

    ray.setOrigin(chemkit::Point3f(1, 2, 3));
    QCOMPARE(ray.origin(), chemkit::Point3f(1, 2, 3));
}

void GraphicsRayTest::setDirection()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.direction(), chemkit::Vector3f(0, 0, -1));

    ray.setDirection(chemkit::Vector3f(1, 0, 0));
    QCOMPARE(ray.direction(), chemkit::Vector3f(1, 0, 0));

    ray.setDirection(chemkit::Vector3f(0, 4, 0));
    QCOMPARE(ray.direction(), chemkit::Vector3f(0, 1, 0));
}

void GraphicsRayTest::intersectsSphere()
{
    chemkit::GraphicsRay ray(chemkit::Point3f(0, 0, 0), chemkit::Vector3f(0, 1, 0));
    float distance;

    QCOMPARE(ray.intersectsSphere(chemkit::Point3f(0, 2, 0), 1.0), true);

    QCOMPARE(ray.intersectsSphere(chemkit::Point3f(0, 2, 0), 1.0, &distance), true);
    QCOMPARE(distance, 1.0f);

    QCOMPARE(ray.intersectsSphere(chemkit::Point3f(0, 0, 0), 2.0, &distance), true);
    QCOMPARE(distance, 2.0f);

    QCOMPARE(ray.intersectsSphere(chemkit::Point3f(4, 0, 0), 1.5), false);

    ray = chemkit::GraphicsRay(chemkit::Point3f(5, 6, 7), chemkit::Vector3f(-1, 0, 0));

    QCOMPARE(ray.intersectsSphere(chemkit::Point3f(-3, 6, 7), 1.0, &distance), true);
    QCOMPARE(distance, 7.0f);
}

void GraphicsRayTest::intersectsCylinder()
{
    chemkit::GraphicsRay ray(chemkit::Point3f(0, 0, 0), chemkit::Vector3f(0, 1, 0));
    float distance;

    QCOMPARE(ray.intersectsCylinder(chemkit::Point3f(2, 2, 0), chemkit::Point3f(-2, 2, 0), 1.0, &distance), true);
    QCOMPARE(distance, 1.0f);
}

void GraphicsRayTest::pointAt()
{
    chemkit::GraphicsRay ray(chemkit::Point3f(1, 0, 0), chemkit::Vector3f(1, 0, 0));
    QCOMPARE(ray.pointAt(0), ray.origin());
    QCOMPARE(ray.pointAt(1.0), chemkit::Point3f(2, 0, 0));
    QCOMPARE(ray.pointAt(-4.0), chemkit::Point3f(-3, 0, 0));

    ray.setDirection(chemkit::Vector3f(0, 1, 0));
    QCOMPARE(ray.pointAt(2.0), chemkit::Point3f(1, 2, 0));
}

QTEST_APPLESS_MAIN(GraphicsRayTest)
