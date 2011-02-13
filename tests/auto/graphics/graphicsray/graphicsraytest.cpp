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

#include <QtTest>

#include <chemkit/graphicsray.h>

class GraphicsRayTest : public QObject
{
    Q_OBJECT

    private slots:
        void basic();
        void setOrigin();
        void setDirection();
        void intersectsSphere();
        void intersectsCylinder();
        void pointAt();
};

void GraphicsRayTest::basic()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.origin(), chemkit::GraphicsPoint(0, 0, 0));
    QCOMPARE(ray.direction(), chemkit::GraphicsVector(0, 0, -1));

    ray = chemkit::GraphicsRay(chemkit::GraphicsPoint(0, 0, 0), chemkit::GraphicsVector(0, 1, 0));
    QCOMPARE(ray.origin(), chemkit::GraphicsPoint(0, 0, 0));
    QCOMPARE(ray.direction(), chemkit::GraphicsVector(0, 1, 0));
}

void GraphicsRayTest::setOrigin()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.origin(), chemkit::GraphicsPoint(0, 0, 0));

    ray.setOrigin(chemkit::GraphicsPoint(1, 2, 3));
    QCOMPARE(ray.origin(), chemkit::GraphicsPoint(1, 2, 3));
}

void GraphicsRayTest::setDirection()
{
    chemkit::GraphicsRay ray;
    QCOMPARE(ray.direction(), chemkit::GraphicsVector(0, 0, -1));

    ray.setDirection(chemkit::GraphicsVector(1, 0, 0));
    QCOMPARE(ray.direction(), chemkit::GraphicsVector(1, 0, 0));

    ray.setDirection(chemkit::GraphicsVector(0, 4, 0));
    QCOMPARE(ray.direction(), chemkit::GraphicsVector(0, 1, 0));
}

void GraphicsRayTest::intersectsSphere()
{
    chemkit::GraphicsRay ray(chemkit::GraphicsPoint(0, 0, 0), chemkit::GraphicsVector(0, 1, 0));
    chemkit::GraphicsFloat distance;

    QCOMPARE(ray.intersectsSphere(chemkit::GraphicsPoint(0, 2, 0), 1.0), true);

    QCOMPARE(ray.intersectsSphere(chemkit::GraphicsPoint(0, 2, 0), 1.0, &distance), true);
    QCOMPARE(distance, 1.0f);

    QCOMPARE(ray.intersectsSphere(chemkit::GraphicsPoint(0, 0, 0), 2.0, &distance), true);
    QCOMPARE(distance, 2.0f);

    QCOMPARE(ray.intersectsSphere(chemkit::GraphicsPoint(4, 0, 0), 1.5), false);

    ray = chemkit::GraphicsRay(chemkit::GraphicsPoint(5, 6, 7), chemkit::GraphicsVector(-1, 0, 0));

    QCOMPARE(ray.intersectsSphere(chemkit::GraphicsPoint(-3, 6, 7), 1.0, &distance), true);
    QCOMPARE(distance, 7.0f);
}

void GraphicsRayTest::intersectsCylinder()
{
    chemkit::GraphicsRay ray(chemkit::GraphicsPoint(0, 0, 0), chemkit::GraphicsVector(0, 1, 0));
    chemkit::GraphicsFloat distance;

    QCOMPARE(ray.intersectsCylinder(chemkit::GraphicsPoint(2, 2, 0), chemkit::GraphicsPoint(-2, 2, 0), 1.0, &distance), true);
    QCOMPARE(distance, 1.0f);
}

void GraphicsRayTest::pointAt()
{
    chemkit::GraphicsRay ray(chemkit::GraphicsPoint(1, 0, 0), chemkit::GraphicsVector(1, 0, 0));
    QCOMPARE(ray.pointAt(0), ray.origin());
    QCOMPARE(ray.pointAt(1.0), chemkit::GraphicsPoint(2, 0, 0));
    QCOMPARE(ray.pointAt(-4.0), chemkit::GraphicsPoint(-3, 0, 0));

    ray.setDirection(chemkit::GraphicsVector(0, 1, 0));
    QCOMPARE(ray.pointAt(2.0), chemkit::GraphicsPoint(1, 2, 0));
}

QTEST_APPLESS_MAIN(GraphicsRayTest)
#include "graphicsraytest.moc"
