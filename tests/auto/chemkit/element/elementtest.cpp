/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "elementtest.h" 

#include <chemkit/element.h>

void ElementTest::symbol()
{
    QCOMPARE(chemkit::Element().symbol(), std::string());
    QCOMPARE(chemkit::Element(6).symbol(), std::string("C"));
    QCOMPARE(chemkit::Element(1).symbol(), std::string("H"));
    QCOMPARE(chemkit::Element(200).symbol(), std::string());
    QCOMPARE(chemkit::Element(-1).symbol(), std::string());
    QCOMPARE(chemkit::Element(109).symbol(), std::string("Mt"));
}

void ElementTest::name()
{
    QCOMPARE(chemkit::Element().name(), std::string());
    QCOMPARE(chemkit::Element(6).name(), std::string("Carbon"));
    QCOMPARE(chemkit::Element(1).name(), std::string("Hydrogen"));
    QCOMPARE(chemkit::Element(200).name(), std::string());
    QCOMPARE(chemkit::Element(-1).name(), std::string());
    QCOMPARE(chemkit::Element(109).name(), std::string("Meitnerium"));
}

void ElementTest::atomicNumber()
{
    QCOMPARE(chemkit::Element("C").atomicNumber(), chemkit::Element::AtomicNumberType(6));
    QCOMPARE(chemkit::Element("H").atomicNumber(), chemkit::Element::AtomicNumberType(1));
    QCOMPARE(chemkit::Element("").atomicNumber(), chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element(std::string()).atomicNumber(), chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element("X").atomicNumber(), chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element("Xa").atomicNumber(), chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element("Xaa").atomicNumber(), chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element("He").atomicNumber(), chemkit::Element::AtomicNumberType(2));
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

void ElementTest::fromName()
{
    QCOMPARE(chemkit::Element::fromName("Hydrogen").atomicNumber(),
             chemkit::Element::AtomicNumberType(1));
    QCOMPARE(chemkit::Element::fromName("Carbon").atomicNumber(),
             chemkit::Element::AtomicNumberType(6));
    QCOMPARE(chemkit::Element::fromName("InvalidName").atomicNumber(),
             chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element::fromName("").atomicNumber(),
             chemkit::Element::AtomicNumberType(0));
}

void ElementTest::fromSymbol()
{
    QCOMPARE(chemkit::Element::fromSymbol("H").atomicNumber(),
             chemkit::Element::AtomicNumberType(1));
    QCOMPARE(chemkit::Element::fromSymbol("C").atomicNumber(),
             chemkit::Element::AtomicNumberType(6));
    QCOMPARE(chemkit::Element::fromSymbol("Invalid").atomicNumber(),
             chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element::fromSymbol("").atomicNumber(),
             chemkit::Element::AtomicNumberType(0));
    QCOMPARE(chemkit::Element::fromSymbol('N').atomicNumber(),
             chemkit::Element::AtomicNumberType(7));
    QCOMPARE(chemkit::Element::fromSymbol("FeAtom", 2).atomicNumber(),
             chemkit::Element::AtomicNumberType(26));
    QCOMPARE(chemkit::Element::fromSymbol("S", 1).atomicNumber(),
             chemkit::Element::AtomicNumberType(16));
}

void ElementTest::isValidAtomicNumber()
{
    QCOMPARE(chemkit::Element::isValidAtomicNumber(1), true);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(6), true);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(0), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(200), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(-1), false);
    QCOMPARE(chemkit::Element::isValidAtomicNumber(109), true);
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
