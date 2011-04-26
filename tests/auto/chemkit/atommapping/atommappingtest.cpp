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
