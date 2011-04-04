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
