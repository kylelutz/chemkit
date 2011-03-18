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

#include <chemkit/point3.h>

class Point3Test : public QObject
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

void Point3Test::basic()
{
    chemkit::Point3 point;
    QCOMPARE(point.x(), chemkit::Float(0.0));
    QCOMPARE(point.y(), chemkit::Float(0.0));
    QCOMPARE(point.z(), chemkit::Float(0.0));

    point = chemkit::Point3(1, 2, 3);
    QCOMPARE(point.x(), chemkit::Float(1.0));
    QCOMPARE(point.y(), chemkit::Float(2.0));
    QCOMPARE(point.z(), chemkit::Float(3.0));
}

void Point3Test::moveBy()
{
}

void Point3Test::movedBy()
{
}

void Point3Test::isNull()
{
    chemkit::Point3 point;
    QCOMPARE(point.isNull(), true);

    point = chemkit::Point3(1, 2, 3);
    QCOMPARE(point.isNull(), false);

    point = chemkit::Point3(0, 0, 0);
    QCOMPARE(point.isNull(), true);
}

void Point3Test::distance()
{
    chemkit::Point3 a, b;
    QCOMPARE(a.distance(b), chemkit::Float(0.0));

    a = chemkit::Point3(1, 0, 0);
    b = chemkit::Point3(3, 0, 0);
    QCOMPARE(a.distance(b), chemkit::Float(2.0));
}

void Point3Test::midpoint()
{
    chemkit::Point3 a(1, 0, 0);
    chemkit::Point3 b(3, 0, 0);

    chemkit::Point3 midpoint = a.midpoint(b);
    QCOMPARE(midpoint, chemkit::Point3(2, 0, 0));
}

QTEST_APPLESS_MAIN(Point3Test)
#include "point3test.moc"
