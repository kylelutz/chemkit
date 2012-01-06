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

#include "substructurequerytest.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/substructurequery.h>

void SubstructureQueryTest::molecule()
{
    chemkit::SubstructureQuery query;
    QVERIFY(query.molecule() == 0);

    chemkit::Molecule molecule;
    query.setMolecule(&molecule);
    QVERIFY(query.molecule() == &molecule);

    query.setMolecule(0);
    QVERIFY(query.molecule() == 0);
}

void SubstructureQueryTest::mapping()
{
    chemkit::Molecule methanol;
    chemkit::Atom *methanol_C1 = methanol.addAtom("C");
    chemkit::Atom *methanol_O2 = methanol.addAtom("O");
    chemkit::Atom *methanol_H3 = methanol.addAtom("H");
    methanol.addBond(methanol_C1, methanol_O2);
    methanol.addBond(methanol_O2, methanol_H3);

    chemkit::Molecule ethanol;
    chemkit::Atom *ethanol_C1 = ethanol.addAtom("C");
    chemkit::Atom *ethanol_C2 = ethanol.addAtom("C");
    chemkit::Atom *ethanol_O3 = ethanol.addAtom("O");
    chemkit::Atom *ethanol_H4 = ethanol.addAtom("H");
    ethanol.addBond(ethanol_C1, ethanol_C2);
    ethanol.addBond(ethanol_C2, ethanol_O3);
    ethanol.addBond(ethanol_O3, ethanol_H4);

    chemkit::SubstructureQuery query(&methanol);
    std::map<chemkit::Atom *, chemkit::Atom *> mapping = query.mapping(&ethanol);
    QVERIFY(mapping.size() == 2);
    QVERIFY(mapping[methanol_C1] == ethanol_C2);
    QVERIFY(mapping[methanol_O2] == ethanol_O3);

    query.setFlags(chemkit::SubstructureQuery::CompareHydrogens);
    mapping = query.mapping(&ethanol);
    QVERIFY(mapping.size() == 3);
    QVERIFY(mapping[methanol_C1] == ethanol_C2);
    QVERIFY(mapping[methanol_O2] == ethanol_O3);
    QVERIFY(mapping[methanol_H3] == ethanol_H4);
}

void SubstructureQueryTest::matches()
{
    chemkit::SubstructureQuery query;

    chemkit::Molecule empty1;
    chemkit::Molecule empty2;
    query.setMolecule(&empty1);
    QCOMPARE(query.matches(&empty2), true);
    query.setMolecule(&empty2);
    QCOMPARE(query.matches(&empty1), true);

    chemkit::Molecule methane("C", "smiles");
    chemkit::Molecule ethane("CC", "smiles");
    chemkit::Molecule propane("CCC", "smiles");
    chemkit::Molecule benzene("c1ccccc1", "smiles");
    chemkit::Molecule phenol("c1ccccc1O", "smiles");

    query.setMolecule(&methane);
    QCOMPARE(query.matches(&methane), true);
    QCOMPARE(query.matches(&ethane), true);
    QCOMPARE(query.matches(&propane), true);
    QCOMPARE(query.matches(&benzene), true);
    QCOMPARE(query.matches(&phenol), true);

    query.setMolecule(&ethane);
    QCOMPARE(query.matches(&methane), false);
    QCOMPARE(query.matches(&ethane), true);
    QCOMPARE(query.matches(&propane), true);
    QCOMPARE(query.matches(&benzene), true);
    QCOMPARE(query.matches(&phenol), true);

    query.setMolecule(&propane);
    QCOMPARE(query.matches(&methane), false);
    QCOMPARE(query.matches(&ethane), false);
    QCOMPARE(query.matches(&propane), true);
    QCOMPARE(query.matches(&benzene), false);
    QCOMPARE(query.matches(&phenol), false);

    query.setMolecule(&benzene);
    QCOMPARE(query.matches(&methane), false);
    QCOMPARE(query.matches(&ethane), false);
    QCOMPARE(query.matches(&propane), false);
    QCOMPARE(query.matches(&benzene), true);
    QCOMPARE(query.matches(&phenol), true);

    query.setMolecule(&phenol);
    QCOMPARE(query.matches(&methane), false);
    QCOMPARE(query.matches(&ethane), false);
    QCOMPARE(query.matches(&propane), false);
    QCOMPARE(query.matches(&benzene), false);
    QCOMPARE(query.matches(&phenol), true);
}

void SubstructureQueryTest::find()
{
    chemkit::Molecule alanine;
    chemkit::Atom *alanine_C1 = alanine.addAtom("C");
    chemkit::Atom *alanine_C2 = alanine.addAtom("C");
    chemkit::Atom *alanine_C3 = alanine.addAtom("C");
    chemkit::Atom *alanine_N4 = alanine.addAtom("N");
    chemkit::Atom *alanine_O5 = alanine.addAtom("O");
    chemkit::Atom *alanine_O6 = alanine.addAtom("O");
    alanine.addBond(alanine_C1, alanine_C2);
    alanine.addBond(alanine_C1, alanine_C3);
    alanine.addBond(alanine_C1, alanine_N4);
    alanine.addBond(alanine_C2, alanine_O5, chemkit::Bond::Double);
    alanine.addBond(alanine_C2, alanine_O6);

    chemkit::Molecule carbonyl;
    chemkit::Atom *carbonyl_C1 = carbonyl.addAtom("C");
    chemkit::Atom *carbonyl_O2 = carbonyl.addAtom("O");
    carbonyl.addBond(carbonyl_C1, carbonyl_O2, chemkit::Bond::Double);

    chemkit::Molecule carboxyl;
    chemkit::Atom *carboxyl_C1 = carboxyl.addAtom("C");
    chemkit::Atom *carboxyl_O2 = carboxyl.addAtom("O");
    chemkit::Atom *carboxyl_O3 = carboxyl.addAtom("O");
    carboxyl.addBond(carboxyl_C1, carboxyl_O2, chemkit::Bond::Double);
    carboxyl.addBond(carboxyl_C1, carboxyl_O3);

    chemkit::SubstructureQuery query(&carboxyl);
    chemkit::Moiety carboxylMoiety = query.find(&alanine);
    QVERIFY(carboxylMoiety.molecule() == &alanine);
    QCOMPARE(carboxylMoiety.atomCount(), size_t(3));
    QVERIFY(carboxylMoiety.atom(0) == alanine_C2);
    QVERIFY(carboxylMoiety.atom(1) == alanine_O5);
    QVERIFY(carboxylMoiety.atom(2) == alanine_O6);

    query.setMolecule(&carbonyl);
    chemkit::Moiety carbonylMoiety = query.find(&carboxyl);
    QVERIFY(carbonylMoiety.molecule() == &carboxyl);
    QCOMPARE(carbonylMoiety.atomCount(), size_t(2));
    QVERIFY(carbonylMoiety.atom(0) == carboxyl_C1);
    QVERIFY(carbonylMoiety.atom(1) == carboxyl_O2);

    query.setMolecule(&carboxyl);
    carboxylMoiety = query.find(&carbonyl);
    QCOMPARE(carboxylMoiety.isEmpty(), true);
}

QTEST_APPLESS_MAIN(SubstructureQueryTest)
