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

#include "formulatest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void FormulaTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::LineFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "formula") != formats.end());
}

void FormulaTest::read()
{
    chemkit::LineFormat *formulaFormat = chemkit::LineFormat::create("formula");
    QVERIFY(formulaFormat != 0);

    // empty
    chemkit::Molecule *empty = formulaFormat->read("");
    QCOMPARE(empty->isEmpty(), true);
    delete empty;

    // hexane
    chemkit::Molecule *hexane = formulaFormat->read("C6H14");
    QCOMPARE(hexane->atomCount(), 20);
    QCOMPARE(hexane->atomCount(chemkit::Atom::Carbon), 6);
    QCOMPARE(hexane->atomCount(chemkit::Atom::Hydrogen), 14);
    delete hexane;

    // water
    chemkit::Molecule *water = formulaFormat->read("H2O");
    QCOMPARE(water->atomCount(), 3);
    QCOMPARE(water->atomCount(chemkit::Atom::Hydrogen), 2);
    QCOMPARE(water->atomCount(chemkit::Atom::Oxygen), 1);
    delete water;

    // atp
    chemkit::Molecule *atp = formulaFormat->read("C10H16N5O13P3");
    QCOMPARE(atp->atomCount(), 47);
    QCOMPARE(atp->atomCount(chemkit::Atom::Carbon), 10);
    QCOMPARE(atp->atomCount(chemkit::Atom::Hydrogen), 16);
    QCOMPARE(atp->atomCount(chemkit::Atom::Nitrogen), 5);
    QCOMPARE(atp->atomCount(chemkit::Atom::Oxygen), 13);
    QCOMPARE(atp->atomCount(chemkit::Atom::Phosphorus), 3);
    delete atp;

    delete formulaFormat;
}

void FormulaTest::write()
{
    chemkit::LineFormat *formulaFormat = chemkit::LineFormat::create("formula");
    QVERIFY(formulaFormat != 0);

    // empty
    chemkit::Molecule empty;
    QCOMPARE(formulaFormat->write(&empty), std::string());

    // water
    chemkit::Molecule water;
    water.addAtom("H");
    water.addAtom("H");
    water.addAtom("O");
    QCOMPARE(formulaFormat->write(&water), std::string("H2O"));

    delete formulaFormat;
}

QTEST_APPLESS_MAIN(FormulaTest)
