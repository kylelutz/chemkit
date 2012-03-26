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

#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/substructurequery.h>

void SubstructureQueryTest::molecule()
{
    chemkit::SubstructureQuery query;
    QVERIFY(query.molecule() == 0);

    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    query.setMolecule(molecule);
    QVERIFY(query.molecule() == molecule);

    query.setMolecule(boost::shared_ptr<chemkit::Molecule>());
    QVERIFY(query.molecule() == boost::shared_ptr<chemkit::Molecule>());
}

void SubstructureQueryTest::mapping()
{
    boost::shared_ptr<chemkit::Molecule> methanol(new chemkit::Molecule);
    chemkit::Atom *methanol_C1 = methanol->addAtom("C");
    chemkit::Atom *methanol_O2 = methanol->addAtom("O");
    chemkit::Atom *methanol_H3 = methanol->addAtom("H");
    methanol->addBond(methanol_C1, methanol_O2);
    methanol->addBond(methanol_O2, methanol_H3);

    boost::shared_ptr<chemkit::Molecule> ethanol(new chemkit::Molecule);
    chemkit::Atom *ethanol_C1 = ethanol->addAtom("C");
    chemkit::Atom *ethanol_C2 = ethanol->addAtom("C");
    chemkit::Atom *ethanol_O3 = ethanol->addAtom("O");
    chemkit::Atom *ethanol_H4 = ethanol->addAtom("H");
    ethanol->addBond(ethanol_C1, ethanol_C2);
    ethanol->addBond(ethanol_C2, ethanol_O3);
    ethanol->addBond(ethanol_O3, ethanol_H4);

    chemkit::SubstructureQuery query(methanol);
    std::map<chemkit::Atom *, chemkit::Atom *> mapping = query.mapping(ethanol.get());
    QVERIFY(mapping.size() == 2);
    QVERIFY(mapping[methanol_C1] == ethanol_C2);
    QVERIFY(mapping[methanol_O2] == ethanol_O3);

    query.setFlags(chemkit::SubstructureQuery::CompareHydrogens);
    mapping = query.mapping(ethanol.get());
    QVERIFY(mapping.size() == 3);
    QVERIFY(mapping[methanol_C1] == ethanol_C2);
    QVERIFY(mapping[methanol_O2] == ethanol_O3);
    QVERIFY(mapping[methanol_H3] == ethanol_H4);
}

void SubstructureQueryTest::matches()
{
    chemkit::SubstructureQuery query;

    boost::shared_ptr<chemkit::Molecule> empty1(new chemkit::Molecule);
    boost::shared_ptr<chemkit::Molecule> empty2(new chemkit::Molecule);
    query.setMolecule(empty1);
    QCOMPARE(query.matches(empty2.get()), true);
    query.setMolecule(empty2);
    QCOMPARE(query.matches(empty1.get()), true);

    boost::shared_ptr<chemkit::Molecule> methane = boost::make_shared<chemkit::Molecule>("C", "smiles");
    boost::shared_ptr<chemkit::Molecule> ethane = boost::make_shared<chemkit::Molecule>("CC", "smiles");
    boost::shared_ptr<chemkit::Molecule> propane = boost::make_shared<chemkit::Molecule>("CCC", "smiles");
    boost::shared_ptr<chemkit::Molecule> benzene = boost::make_shared<chemkit::Molecule>("c1ccccc1", "smiles");
    boost::shared_ptr<chemkit::Molecule> phenol = boost::make_shared<chemkit::Molecule>("c1ccccc1O", "smiles");

    query.setMolecule(methane);
    QCOMPARE(query.matches(methane.get()), true);
    QCOMPARE(query.matches(ethane.get()), true);
    QCOMPARE(query.matches(propane.get()), true);
    QCOMPARE(query.matches(benzene.get()), true);
    QCOMPARE(query.matches(phenol.get()), true);

    query.setMolecule(ethane);
    QCOMPARE(query.matches(methane.get()), false);
    QCOMPARE(query.matches(ethane.get()), true);
    QCOMPARE(query.matches(propane.get()), true);
    QCOMPARE(query.matches(benzene.get()), true);
    QCOMPARE(query.matches(phenol.get()), true);

    query.setMolecule(propane);
    QCOMPARE(query.matches(methane.get()), false);
    QCOMPARE(query.matches(ethane.get()), false);
    QCOMPARE(query.matches(propane.get()), true);
    QCOMPARE(query.matches(benzene.get()), false);
    QCOMPARE(query.matches(phenol.get()), false);

    query.setMolecule(benzene);
    QCOMPARE(query.matches(methane.get()), false);
    QCOMPARE(query.matches(ethane.get()), false);
    QCOMPARE(query.matches(propane.get()), false);
    QCOMPARE(query.matches(benzene.get()), true);
    QCOMPARE(query.matches(phenol.get()), true);

    query.setMolecule(phenol);
    QCOMPARE(query.matches(methane.get()), false);
    QCOMPARE(query.matches(ethane.get()), false);
    QCOMPARE(query.matches(propane.get()), false);
    QCOMPARE(query.matches(benzene.get()), false);
    QCOMPARE(query.matches(phenol.get()), true);
}

void SubstructureQueryTest::find()
{
    boost::shared_ptr<chemkit::Molecule> alanine(new chemkit::Molecule);
    chemkit::Atom *alanine_C1 = alanine->addAtom("C");
    chemkit::Atom *alanine_C2 = alanine->addAtom("C");
    chemkit::Atom *alanine_C3 = alanine->addAtom("C");
    chemkit::Atom *alanine_N4 = alanine->addAtom("N");
    chemkit::Atom *alanine_O5 = alanine->addAtom("O");
    chemkit::Atom *alanine_O6 = alanine->addAtom("O");
    alanine->addBond(alanine_C1, alanine_C2);
    alanine->addBond(alanine_C1, alanine_C3);
    alanine->addBond(alanine_C1, alanine_N4);
    alanine->addBond(alanine_C2, alanine_O5, chemkit::Bond::Double);
    alanine->addBond(alanine_C2, alanine_O6);

    boost::shared_ptr<chemkit::Molecule> carbonyl(new chemkit::Molecule);
    chemkit::Atom *carbonyl_C1 = carbonyl->addAtom("C");
    chemkit::Atom *carbonyl_O2 = carbonyl->addAtom("O");
    carbonyl->addBond(carbonyl_C1, carbonyl_O2, chemkit::Bond::Double);

    boost::shared_ptr<chemkit::Molecule> carboxyl(new chemkit::Molecule);
    chemkit::Atom *carboxyl_C1 = carboxyl->addAtom("C");
    chemkit::Atom *carboxyl_O2 = carboxyl->addAtom("O");
    chemkit::Atom *carboxyl_O3 = carboxyl->addAtom("O");
    carboxyl->addBond(carboxyl_C1, carboxyl_O2, chemkit::Bond::Double);
    carboxyl->addBond(carboxyl_C1, carboxyl_O3);

    chemkit::SubstructureQuery query(carboxyl);
    chemkit::Moiety carboxylMoiety = query.find(alanine.get());
    QVERIFY(carboxylMoiety.molecule() == alanine.get());
    QCOMPARE(carboxylMoiety.atomCount(), size_t(3));
    QVERIFY(carboxylMoiety.atom(0) == alanine_C2);
    QVERIFY(carboxylMoiety.atom(1) == alanine_O5);
    QVERIFY(carboxylMoiety.atom(2) == alanine_O6);

    query.setMolecule(carbonyl);
    chemkit::Moiety carbonylMoiety = query.find(carboxyl.get());
    QVERIFY(carbonylMoiety.molecule() == carboxyl.get());
    QCOMPARE(carbonylMoiety.atomCount(), size_t(2));
    QVERIFY(carbonylMoiety.atom(0) == carboxyl_C1);
    QVERIFY(carbonylMoiety.atom(1) == carboxyl_O2);

    query.setMolecule(carboxyl);
    carboxylMoiety = query.find(carbonyl.get());
    QCOMPARE(carboxylMoiety.isEmpty(), true);
}

QTEST_APPLESS_MAIN(SubstructureQueryTest)
