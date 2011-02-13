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

#include <chemkit/matrix.h>

class MatrixTest : public QObject
{
    Q_OBJECT

    private slots:
        void multiply();
        void identity();
};

void MatrixTest::multiply()
{
    chemkit::Matrix a(2, 4);
    a = 2, 4, 8, 9,
        -10, 2, 6, 15;

    chemkit::Matrix b(4, 2);
    b = 8, 9,
        7, 4.5,
        9, 10,
        -2, 4;

    chemkit::Matrix c = a.multiply(b);
    QCOMPARE(c.rowCount(), 2);
    QCOMPARE(c.columnCount(), 2);

    chemkit::Matrix d(2, 2);
    d = 98, 152,
        -42, 39;

    QVERIFY(c == d);
}

void MatrixTest::identity()
{
    chemkit::Matrix identity = chemkit::Matrix::identity(4, 4);
    QCOMPARE(qRound(identity(0, 0)), 1);
    QCOMPARE(qRound(identity(0, 1)), 0);
    QCOMPARE(qRound(identity(2, 2)), 1);
    QCOMPARE(qRound(identity(2, 3)), 0);
}

QTEST_APPLESS_MAIN(MatrixTest)
#include "matrixtest.moc"
