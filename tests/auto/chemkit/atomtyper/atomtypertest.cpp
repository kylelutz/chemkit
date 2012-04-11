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
    QCOMPARE(typer->type(molecule.atom(0)), std::string("C"));
    QCOMPARE(typer->type(molecule.atom(1)), std::string("O"));
    QCOMPARE(typer->type(molecule.atom(2)), std::string("S"));
    QCOMPARE(typer->type(molecule.atom(3)), std::string("Pb"));

    delete typer;
}

void AtomTyperTest::cleanupTestCase()
{
    delete m_plugin;

    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "mock") == typers.end());
}

QTEST_APPLESS_MAIN(AtomTyperTest)
