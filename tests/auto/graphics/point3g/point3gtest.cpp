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

#include <chemkit/point3g.h>

class Point3gTest : public QObject
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

void Point3gTest::xyz()
{
    chemkit::Point3g nullpoint;
    QCOMPARE(nullpoint.x(), 0.0f);
    QCOMPARE(nullpoint.y(), 0.0f);
    QCOMPARE(nullpoint.z(), 0.0f);

    chemkit::Point3g point(-1.5f, 8.0f, 4.1f);
    QCOMPARE(point.x(), -1.5f);
    QCOMPARE(point.y(), 8.0f);
    QCOMPARE(point.z(), 4.1f);
}

void Point3gTest::moveBy()
{
    chemkit::Point3g point(0, 0, 0);
    point.moveBy(2, 0, 0);
    QCOMPARE(point.x(), 2.0f);
    QCOMPARE(point.y(), 0.0f);
    QCOMPARE(point.z(), 0.0f);

    point.moveBy(-3, 8, 9);
    QCOMPARE(point.x(), -1.0f);
    QCOMPARE(point.y(), 8.0f);
    QCOMPARE(point.z(), 9.0f);
}

void Point3gTest::movedBy()
{
    chemkit::Point3g point(-1, 0, 4);
    chemkit::Point3g movedPoint = point.movedBy(0, 5, -1);

    // check that point is the same
    QCOMPARE(point.x(), -1.0f);
    QCOMPARE(point.y(), 0.0f);
    QCOMPARE(point.z(), 4.0f);

    // check that the movedPoint moved
    QCOMPARE(movedPoint.x(), -1.0f);
    QCOMPARE(movedPoint.y(), 5.0f);
    QCOMPARE(movedPoint.z(), 3.0f);
}

void Point3gTest::isNull()
{
    chemkit::Point3g point;
    QCOMPARE(point.isNull(), true);

    point = chemkit::Point3g(1, 2, 3);
    QCOMPARE(point.isNull(), false);

    point.moveBy(-1, -2, -3);
    QCOMPARE(point.isNull(), true);
}

void Point3gTest::distance()
{
    chemkit::Point3g a(0, 0, 0);
    chemkit::Point3g b(0, 0, 0);
    QCOMPARE(a.distance(b), 0.0f);

    a = chemkit::Point3g(2, 0, 0);
    QCOMPARE(a.distance(b), 2.0f);
}

void Point3gTest::midpoint()
{
    chemkit::Point3g a(0, 0, 0);
    chemkit::Point3g b(0, 0, 0);
    QVERIFY(a.midpoint(b) == chemkit::Point3g(0, 0, 0));

    a = chemkit::Point3g(0, 4, 0);
    QVERIFY(a.midpoint(b) == chemkit::Point3g(0, 2, 0));

    b = chemkit::Point3g(0, 2, 0);
    QVERIFY(b.midpoint(a) == chemkit::Point3g(0, 3, 0));
}

QTEST_APPLESS_MAIN(Point3gTest)
#include "point3gtest.moc"
