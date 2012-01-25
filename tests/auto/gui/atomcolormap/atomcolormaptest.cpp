/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "atomcolormaptest.h"

#include <chemkit/atomcolormap.h>

void AtomColorMapTest::color()
{
    chemkit::AtomColorMap map;
    QCOMPARE(map.color("C"), QColor());

    map.setColor("C", Qt::red);
    QCOMPARE(map.color("C"), QColor(Qt::red));
}

void AtomColorMapTest::defaultColor()
{
    chemkit::AtomColorMap colorMap;
    QCOMPARE(colorMap.defaultColor(), QColor());

    colorMap.setDefaultColor(Qt::blue);
    QCOMPARE(colorMap.defaultColor(), QColor(Qt::blue));

    chemkit::AtomColorMap defaultColorMap(chemkit::AtomColorMap::DefaultColorScheme);
    QCOMPARE(defaultColorMap.defaultColor(), QColor(255, 20, 147));

    chemkit::AtomColorMap rasmolColorMap(chemkit::AtomColorMap::RasmolColorScheme);
    QCOMPARE(rasmolColorMap.defaultColor(), QColor(255, 20, 147));

    chemkit::AtomColorMap pymolColorMap(chemkit::AtomColorMap::PymolColorScheme);
    QCOMPARE(pymolColorMap.defaultColor(), QColor(255, 20, 147));

    chemkit::AtomColorMap jmolColorMap(chemkit::AtomColorMap::JmolColorScheme);
    QCOMPARE(jmolColorMap.defaultColor(), QColor(255, 20, 147));
}

void AtomColorMapTest::carbonColor()
{
    chemkit::AtomColorMap colorMap;
    QCOMPARE(colorMap.color("C"), QColor());

    chemkit::AtomColorMap defaultColorMap(chemkit::AtomColorMap::DefaultColorScheme);
    QCOMPARE(defaultColorMap.color("C"), QColor(80, 80, 80));

    chemkit::AtomColorMap rasmolColorMap(chemkit::AtomColorMap::RasmolColorScheme);
    QCOMPARE(rasmolColorMap.color("C"), QColor(200, 200, 200));

    chemkit::AtomColorMap pymolColorMap(chemkit::AtomColorMap::PymolColorScheme);
    QCOMPARE(pymolColorMap.color("C"), QColor(51, 255, 51));

    chemkit::AtomColorMap jmolColorMap(chemkit::AtomColorMap::JmolColorScheme);
    QCOMPARE(jmolColorMap.color("C"), QColor(144, 144, 144));
}

QTEST_MAIN(AtomColorMapTest)
