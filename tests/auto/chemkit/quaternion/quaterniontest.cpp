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

#include <chemkit/quaternion.h>

class QuaternionTest : public QObject
{
    Q_OBJECT

    private slots:
        void basic();
        void add();
        void subtract();
        void multiply();
        void conjugate();
};

void QuaternionTest::basic()
{
    chemkit::Quaternion q(1, 2, 3, 4);
    QCOMPARE(q.x(), 1.0);
    QCOMPARE(q.y(), 2.0);
    QCOMPARE(q.z(), 3.0);
    QCOMPARE(q.r(), 4.0);

    q = chemkit::Quaternion(5, 6, 7, 8);
    QCOMPARE(q.x(), 5.0);
    QCOMPARE(q.y(), 6.0);
    QCOMPARE(q.z(), 7.0);
    QCOMPARE(q.r(), 8.0);
}

void QuaternionTest::add()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.add(b) == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(a + b == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(b.add(a) == chemkit::Quaternion(6, 8, 10, 12));
    QVERIFY(b + a == chemkit::Quaternion(6, 8, 10, 12));
}

void QuaternionTest::subtract()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.subtract(b) == chemkit::Quaternion(-2, -2, -2, -2));
    QVERIFY(a - b == chemkit::Quaternion(-2, -2, -2, -2));
    QVERIFY(b.subtract(a) == chemkit::Quaternion(2, 2, 2, 2));
    QVERIFY(b - a == chemkit::Quaternion(2, 2, 2, 2));
}

void QuaternionTest::multiply()
{
    chemkit::Quaternion a(2, 3, 4, 5);
    chemkit::Quaternion b(4, 5, 6, 7);
    QVERIFY(a.multiply(b) == chemkit::Quaternion(32, 50, 56, -12));
    QVERIFY(a * b == chemkit::Quaternion(32, 50, 56, -12));
    QVERIFY(b.multiply(a) == chemkit::Quaternion(36, 42, 60, -12));
    QVERIFY(b * a == chemkit::Quaternion(36, 42, 60, -12));
}

void QuaternionTest::conjugate()
{
    chemkit::Quaternion q(1, 2, 3, 4);
    QVERIFY(q.conjugate() == chemkit::Quaternion(-1, -2, -3, 4));

    q = chemkit::Quaternion(-2, -4, -6, -8);
    QVERIFY(q.conjugate() == chemkit::Quaternion(2, 4, 6, -8));
}

QTEST_APPLESS_MAIN(QuaternionTest)
#include "quaterniontest.moc"
