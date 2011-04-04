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

#include "atommappingtest.h"

#include <chemkit/molecule.h>
#include <chemkit/atommapping.h>

void AtomMappingTest::basic()
{
    chemkit::Molecule source;
    chemkit::Molecule target;
    chemkit::AtomMapping mapping(&source, &target);
    QVERIFY(mapping.source() == &source);
    QVERIFY(mapping.target() == &target);

    mapping = chemkit::AtomMapping(&source, &source);
    QVERIFY(mapping.source() == &source);
    QVERIFY(mapping.target() == &source);
}

void AtomMappingTest::size()
{
    chemkit::Molecule source;
    chemkit::Molecule target;
    chemkit::AtomMapping mapping(&source, &target);
    QCOMPARE(mapping.size(), 0);

    chemkit::Atom *source_H1 = source.addAtom("H");
    chemkit::Atom *target_H1 = target.addAtom("H");
    mapping.add(source_H1, target_H1);
    QCOMPARE(mapping.size(), 1);

    mapping.add(source_H1, target_H1);
    QCOMPARE(mapping.size(), 1);

    chemkit::Atom *source_H2 = source.addAtom("H");
    chemkit::Atom *target_H2 = target.addAtom("H");
    mapping.add(source_H2, target_H2);
    QCOMPARE(mapping.size(), 2);

    mapping.clear();
    QCOMPARE(mapping.size(), 0);
}

void AtomMappingTest::isEmpty()
{
    chemkit::Molecule source;
    chemkit::Molecule target;
    chemkit::AtomMapping mapping(&source, &target);
    QCOMPARE(mapping.isEmpty(), true);

    chemkit::Atom *source_H1 = source.addAtom("H");
    chemkit::Atom *target_H1 = target.addAtom("H");
    mapping.add(source_H1, target_H1);
    QCOMPARE(mapping.isEmpty(), false);

    mapping.clear();
    QCOMPARE(mapping.isEmpty(), true);
}

void AtomMappingTest::map()
{
    chemkit::Molecule source;
    chemkit::Molecule target;
    chemkit::AtomMapping mapping(&source, &target);
    chemkit::Atom *source_H1 = source.addAtom("H");
    chemkit::Atom *target_H1 = target.addAtom("H");
    mapping.add(source_H1, target_H1);
    QVERIFY(mapping.map(source_H1) == target_H1);
    QVERIFY(mapping.map(target_H1) == source_H1);

    chemkit::Atom *source_H2 = source.addAtom("H");
    chemkit::Atom *target_H2 = target.addAtom("H");
    QVERIFY(mapping.map(source_H2) == 0);
    QVERIFY(mapping.map(target_H2) == 0);

    mapping.add(source_H2, target_H2);
    QVERIFY(mapping.map(source_H2) == target_H2);
    QVERIFY(mapping.map(target_H2) == source_H2);

    mapping.remove(source_H1);
    QVERIFY(mapping.map(source_H1) == 0);
    QVERIFY(mapping.map(target_H1) == 0);
}

QTEST_APPLESS_MAIN(AtomMappingTest)
