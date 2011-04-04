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

#include "elementtest.h" 

#include <chemkit/element.h>

void ElementTest::symbol()
{
    QCOMPARE(chemkit::Element(6).symbol(), std::string("C"));
    QCOMPARE(chemkit::Element(1).symbol(), std::string("H"));
    QCOMPARE(chemkit::Element(0).symbol(), std::string());
    QCOMPARE(chemkit::Element(500).symbol(), std::string());
    QCOMPARE(chemkit::Element(-1).symbol(), std::string());
    QCOMPARE(chemkit::Element(109).symbol(), std::string("Mt"));
    QCOMPARE(chemkit::Element(110).symbol(), std::string());
}

void ElementTest::name()
{
    QCOMPARE(chemkit::Element(6).name(), std::string("Carbon"));
    QCOMPARE(chemkit::Element(1).name(), std::string("Hydrogen"));
    QCOMPARE(chemkit::Element(0).name(), std::string());
    QCOMPARE(chemkit::Element(500).name(), std::string());
    QCOMPARE(chemkit::Element(-1).name(), std::string());
    QCOMPARE(chemkit::Element(109).name(), std::string("Meitnerium"));
    QCOMPARE(chemkit::Element(110).name(), std::string());
}

void ElementTest::atomicNumber()
{
    QCOMPARE(chemkit::Element("C").atomicNumber(), 6);
    QCOMPARE(chemkit::Element("H").atomicNumber(), 1);
    QCOMPARE(chemkit::Element("").atomicNumber(), 0);
    QCOMPARE(chemkit::Element(std::string()).atomicNumber(), 0);
    QCOMPARE(chemkit::Element("X").atomicNumber(), 0);
    QCOMPARE(chemkit::Element("Xa").atomicNumber(), 0);
    QCOMPARE(chemkit::Element("Xaa").atomicNumber(), 0);
    QCOMPARE(chemkit::Element("He").atomicNumber(), 2);
}

void ElementTest::mass()
{
    QCOMPARE(qRound(chemkit::Element("H").mass()), 1);
    QCOMPARE(qRound(chemkit::Element("C").mass()), 12);
    QCOMPARE(qRound(chemkit::Element("X").mass()), 0);
}

void ElementTest::electronegativity()
{
    QCOMPARE(qRound(chemkit::Element("H").electronegativity()), 2);
    QCOMPARE(qRound(chemkit::Element("C").electronegativity()), 3);
    QCOMPARE(qRound(chemkit::Element("X").electronegativity()), 0);
}

void ElementTest::covalentRadius()
{
    QCOMPARE(qRound(chemkit::Element("K").covalentRadius()), 2);
    QCOMPARE(qRound(chemkit::Element("Cl").covalentRadius()), 1);
    QCOMPARE(qRound(chemkit::Element("X").covalentRadius()), 0);
}

void ElementTest::vanDerWaalsRadius()
{
    QCOMPARE(qRound(chemkit::Element("H").vanDerWaalsRadius()), 1);
    QCOMPARE(qRound(chemkit::Element("C").vanDerWaalsRadius()), 2);
    QCOMPARE(qRound(chemkit::Element("X").vanDerWaalsRadius()), 0);
}

void ElementTest::expectedValence()
{
    QCOMPARE(chemkit::Element(1).expectedValence(), 1);
    QCOMPARE(chemkit::Element(6).expectedValence(), 4);
    QCOMPARE(chemkit::Element(7).expectedValence(), 3);
    QCOMPARE(chemkit::Element(8).expectedValence(), 2);
}

void ElementTest::isMetal()
{
    QCOMPARE(chemkit::Element("H").isMetal(), false);
    QCOMPARE(chemkit::Element("H").isNonmetal(), true);
    QCOMPARE(chemkit::Element("C").isMetal(), false);
    QCOMPARE(chemkit::Element("C").isNonmetal(), true);
    QCOMPARE(chemkit::Element("Li").isMetal(), true);
    QCOMPARE(chemkit::Element("Li").isNonmetal(), false);
}

void ElementTest::isValidAtomicNumber()
{
    QCOMPARE(chemkit::Element::isValidAtomicNumber(1), true);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(6), true);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(0), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(500), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(-1), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(109), true);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(110), false);
}

void ElementTest::isValidSymbol()
{
    QCOMPARE(chemkit::Element::isValidSymbol("C"), true);
    QCOMPARE(chemkit::Element::isValidSymbol("c"), false);
    QCOMPARE(chemkit::Element::isValidSymbol("H"), true);
    QCOMPARE(chemkit::Element::isValidSymbol("h"), false);
    QCOMPARE(chemkit::Element::isValidSymbol("X"), false);
    QCOMPARE(chemkit::Element::isValidSymbol("Xa"), false);
    QCOMPARE(chemkit::Element::isValidSymbol("Xaa"), false);
    QCOMPARE(chemkit::Element::isValidSymbol(""), false);
    QCOMPARE(chemkit::Element::isValidSymbol(std::string()), false);
}

QTEST_APPLESS_MAIN(ElementTest)
