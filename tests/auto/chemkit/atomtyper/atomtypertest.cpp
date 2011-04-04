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

#include "atomtypertest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>

#include "mockatomtyper.h"

void AtomTyperTest::initTestCase()
{
    m_plugin = new MockAtomTyperPlugin;

    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "mock") != typers.end());
}

void AtomTyperTest::create()
{
    chemkit::AtomTyper *typer = chemkit::AtomTyper::create("mock");
    QVERIFY(typer != 0);
    delete typer;

    chemkit::AtomTyper *invalidTyper = chemkit::AtomTyper::create("invalid_name");
    QVERIFY(invalidTyper == 0);
}

void AtomTyperTest::name()
{
    chemkit::AtomTyper *typer = chemkit::AtomTyper::create("mock");
    QVERIFY(typer != 0);
    QCOMPARE(typer->name(), std::string("mock"));
    delete typer;
}

void AtomTyperTest::molecule()
{
    chemkit::AtomTyper *typer = chemkit::AtomTyper::create("mock");
    QVERIFY(typer != 0);
    QVERIFY(typer->molecule() == 0);

    chemkit::Molecule molecule;
    typer->setMolecule(&molecule);
    QVERIFY(typer->molecule() == &molecule);
}

void AtomTyperTest::type()
{
    chemkit::AtomTyper *typer = chemkit::AtomTyper::create("mock");
    QVERIFY(typer != 0);

    chemkit::Molecule molecule;
    molecule.addAtom("C");
    molecule.addAtom("O");
    molecule.addAtom("S");
    molecule.addAtom("Pb");

    typer->setMolecule(&molecule);
    QCOMPARE(typer->type(0).toString(), QString("C"));
    QCOMPARE(typer->type(molecule.atom(1)).toString(), QString("O"));
    QCOMPARE(typer->typeString(2), std::string("S"));
    QCOMPARE(typer->typeString(molecule.atom(3)), std::string("Pb"));

    delete typer;
}

void AtomTyperTest::cleanupTestCase()
{
    delete m_plugin;

    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "mock") == typers.end());
}

QTEST_APPLESS_MAIN(AtomTyperTest)
