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

#include <chemkit/alphashape.h>

class AlphaShapeTest : public QObject
{
    Q_OBJECT

    private slots:
        void alphaValue();
};

void AlphaShapeTest::alphaValue()
{
    std::vector<chemkit::Point3> points;
    chemkit::AlphaShape alphaShape(points);

    // default alpha value is 0
    QCOMPARE(alphaShape.alphaValue(), chemkit::Float(0.0));

    alphaShape.setAlphaValue(1.8);
    QCOMPARE(alphaShape.alphaValue(), chemkit::Float(1.8));
}

QTEST_APPLESS_MAIN(AlphaShapeTest)
#include "alphashapetest.moc"
