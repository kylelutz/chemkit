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

#include <chemkit/point.h>

class PointTest : public QObject
{
    Q_OBJECT

    private slots:
        void basic();
        void moveBy();
        void movedBy();
        void isNull();
        void distance();
        void midpoint();
};

void PointTest::basic()
{
    chemkit::Point point;
    QCOMPARE(point.x(), chemkit::Float(0.0));
    QCOMPARE(point.y(), chemkit::Float(0.0));
    QCOMPARE(point.z(), chemkit::Float(0.0));

    point = chemkit::Point(1, 2, 3);
    QCOMPARE(point.x(), chemkit::Float(1.0));
    QCOMPARE(point.y(), chemkit::Float(2.0));
    QCOMPARE(point.z(), chemkit::Float(3.0));
}

void PointTest::moveBy()
{
}

void PointTest::movedBy()
{
}

void PointTest::isNull()
{
    chemkit::Point point;
    QCOMPARE(point.isNull(), true);

    point = chemkit::Point(1, 2, 3);
    QCOMPARE(point.isNull(), false);

    point = chemkit::Point(0, 0, 0);
    QCOMPARE(point.isNull(), true);
}

void PointTest::distance()
{
    chemkit::Point a, b;
    QCOMPARE(a.distance(b), chemkit::Float(0.0));

    a = chemkit::Point(1, 0, 0);
    b = chemkit::Point(3, 0, 0);
    QCOMPARE(a.distance(b), chemkit::Float(2.0));
}

void PointTest::midpoint()
{
    chemkit::Point a(1, 0, 0);
    chemkit::Point b(3, 0, 0);

    chemkit::Point midpoint = a.midpoint(b);
    QCOMPARE(midpoint, chemkit::Point(2, 0, 0));
}

QTEST_APPLESS_MAIN(PointTest)
#include "pointtest.moc"
