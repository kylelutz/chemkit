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

#include <chemkit/graphicspoint.h>

class GraphicsPointTest : public QObject
{
    Q_OBJECT

    private slots:
        void xyz();
        void moveBy();
        void movedBy();
        void isNull();
        void distance();
        void midpoint();
};

void GraphicsPointTest::xyz()
{
    chemkit::GraphicsPoint nullpoint;
    QCOMPARE(nullpoint.x(), 0.0f);
    QCOMPARE(nullpoint.y(), 0.0f);
    QCOMPARE(nullpoint.z(), 0.0f);

    chemkit::GraphicsPoint point(-1.5f, 8.0f, 4.1f);
    QCOMPARE(point.x(), -1.5f);
    QCOMPARE(point.y(), 8.0f);
    QCOMPARE(point.z(), 4.1f);
}

void GraphicsPointTest::moveBy()
{
    chemkit::GraphicsPoint point(0, 0, 0);
    point.moveBy(2, 0, 0);
    QCOMPARE(point.x(), 2.0f);
    QCOMPARE(point.y(), 0.0f);
    QCOMPARE(point.z(), 0.0f);

    point.moveBy(-3, 8, 9);
    QCOMPARE(point.x(), -1.0f);
    QCOMPARE(point.y(), 8.0f);
    QCOMPARE(point.z(), 9.0f);
}

void GraphicsPointTest::movedBy()
{
    chemkit::GraphicsPoint point(-1, 0, 4);
    chemkit::GraphicsPoint movedPoint = point.movedBy(0, 5, -1);

    // check that point is the same
    QCOMPARE(point.x(), -1.0f);
    QCOMPARE(point.y(), 0.0f);
    QCOMPARE(point.z(), 4.0f);

    // check that the movedPoint moved
    QCOMPARE(movedPoint.x(), -1.0f);
    QCOMPARE(movedPoint.y(), 5.0f);
    QCOMPARE(movedPoint.z(), 3.0f);
}

void GraphicsPointTest::isNull()
{
    chemkit::GraphicsPoint point;
    QCOMPARE(point.isNull(), true);

    point = chemkit::GraphicsPoint(1, 2, 3);
    QCOMPARE(point.isNull(), false);

    point.moveBy(-1, -2, -3);
    QCOMPARE(point.isNull(), true);
}

void GraphicsPointTest::distance()
{
    chemkit::GraphicsPoint a(0, 0, 0);
    chemkit::GraphicsPoint b(0, 0, 0);
    QCOMPARE(a.distance(b), 0.0f);

    a = chemkit::GraphicsPoint(2, 0, 0);
    QCOMPARE(a.distance(b), 2.0f);
}

void GraphicsPointTest::midpoint()
{
    chemkit::GraphicsPoint a(0, 0, 0);
    chemkit::GraphicsPoint b(0, 0, 0);
    QVERIFY(a.midpoint(b) == chemkit::GraphicsPoint(0, 0, 0));

    a = chemkit::GraphicsPoint(0, 4, 0);
    QVERIFY(a.midpoint(b) == chemkit::GraphicsPoint(0, 2, 0));

    b = chemkit::GraphicsPoint(0, 2, 0);
    QVERIFY(b.midpoint(a) == chemkit::GraphicsPoint(0, 3, 0));
}

QTEST_APPLESS_MAIN(GraphicsPointTest)
#include "graphicspointtest.moc"
