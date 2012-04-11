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

#include "moleculetest.h"

#include <boost/bind.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/cartesiancoordinates.h>

void MoleculeTest::name()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.name(), std::string());

    molecule.addAtom("C");
    molecule.addAtom("O");
    QCOMPARE(molecule.name(), std::string());

    molecule.setName("carbonAndOxygen");
    QCOMPARE(molecule.name(), std::string("carbonAndOxygen"));

    molecule.setName(std::string());
    QCOMPARE(molecule.name(), std::string());
}

void MoleculeTest::formula()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.formula(), std::string());

    molecule.addAtom("O");
    molecule.addAtom("O");
    QCOMPARE(molecule.formula(), std::string("O2"));

    molecule.addAtom("Ac");
    QCOMPARE(molecule.formula(), std::string("AcO2"));

    molecule.addAtom("C");
    QCOMPARE(molecule.formula(), std::string("CAcO2"));

    molecule.addAtom("H");
    molecule.addAtom("H");
    QCOMPARE(molecule.formula(), std::string("CH2AcO2"));

    molecule.clear();
    QCOMPARE(molecule.formula(), std::string());
}

void MoleculeTest::mass()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.mass(), chemkit::Real(0.0));

    molecule.addAtom("C");
    QCOMPARE(qRound(molecule.mass()), 12);

    molecule.addAtom("H");
    QCOMPARE(qRound(molecule.mass()), 13);

    molecule.clear();
    QCOMPARE(molecule.mass(), chemkit::Real(0.0));
}

void MoleculeTest::data()
{
    chemkit::Molecule molecule;
    molecule.setData("boilingPoint", 38);
    QCOMPARE(molecule.data("boilingPoint").toInt(), 38);
}

void MoleculeTest::addAtom()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom1;
    chemkit::Atom *atom2;

    atom1 = molecule.addAtom(chemkit::Atom::Carbon);
    QCOMPARE(atom1->atomicNumber(), chemkit::Atom::AtomicNumberType(6));

    atom2 = molecule.addAtom(chemkit::Element());
    QCOMPARE(atom2->element().isValid(), false);

    atom2 = molecule.addAtom(-1);
    QCOMPARE(atom2->element().isValid(), false);

    atom2 = molecule.addAtom(200);
    QCOMPARE(atom2->element().isValid(), false);

    atom1 = molecule.addAtom("C");
    QCOMPARE(atom1->symbol(), std::string("C"));

    atom2 = molecule.addAtom(std::string());
    QCOMPARE(atom2->element().isValid(), false);

    atom2 = molecule.addAtom("");
    QCOMPARE(atom2->element().isValid(), false);

    atom2 = molecule.addAtom("X");
    QCOMPARE(atom2->element().isValid(), false);
}

void MoleculeTest::addAtomCopy()
{
}

void MoleculeTest::removeAtomIf()
{
    chemkit::Molecule ethanol("CCO", "smiles");
    QCOMPARE(ethanol.formula(), std::string("C2H6O"));

    ethanol.removeAtomIf(boost::bind(&chemkit::Atom::isTerminalHydrogen, _1));
    QCOMPARE(ethanol.formula(), std::string("C2O"));

    ethanol.removeAtomIf(boost::bind(&chemkit::Atom::is, _1, chemkit::Atom::Carbon));
    QCOMPARE(ethanol.formula(), std::string("O"));
}

void MoleculeTest::atom()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    QVERIFY(molecule.atom(0) == C1);
    QVERIFY(molecule.atom(1) == C2);

    // check the operator[] method
    QVERIFY(molecule[0] == C1);
    QVERIFY(molecule[1] == C2);
}

void MoleculeTest::addBond()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    QCOMPARE(molecule.bondCount(), size_t(0));

    chemkit::Bond *H1_H2 = molecule.addBond(H1, H2);
    QVERIFY(H1_H2->contains(H1));
    QVERIFY(H1_H2->contains(H2));
    QCOMPARE(H1_H2->order(), chemkit::Bond::BondOrderType(1));
    QCOMPARE(molecule.bondCount(), size_t(1));
    QVERIFY(molecule.bonds()[0] == H1_H2);

    chemkit::Bond *H1_H2_again = molecule.addBond(H1, H2);
    QVERIFY(H1_H2_again == H1_H2);
    QCOMPARE(molecule.bondCount(), size_t(1));

    chemkit::Bond *H1_H1 = molecule.addBond(H1, H1);
    QVERIFY(H1_H1 == 0);
    QCOMPARE(molecule.bondCount(), size_t(1));
}

void MoleculeTest::bond()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    QVERIFY(molecule.bond(C1, C2) == 0);
    QVERIFY(molecule.bond(C2, C3) == 0);
    QVERIFY(molecule.bond(C3, C4) == 0);

    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    QVERIFY(molecule.bond(C1, C2) == C1_C2);
    QVERIFY(molecule.bond(C2, C1) == C1_C2);
    QVERIFY(molecule.bond(C1, C3) == 0);
    QVERIFY(molecule.bond(C3, C2) == 0);

    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    QVERIFY(molecule.bond(C3, C4) == C3_C4);
    QVERIFY(molecule.bond(C4, C3) == C3_C4);
    QVERIFY(molecule.bond(C1, C3) == 0);
    QVERIFY(molecule.bond(C2, C3) == 0);
    QVERIFY(molecule.bond(C2, C4) == 0);

    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    QVERIFY(molecule.bond(C2, C3) == C2_C3);
    QVERIFY(molecule.bond(C3, C2) == C2_C3);
    QVERIFY(molecule.bond(C2, C4) == 0);

    molecule.removeBond(C1_C2);
    QVERIFY(molecule.bond(C1, C2) == 0);
    QVERIFY(molecule.bond(C2, C1) == 0);

    chemkit::Molecule molecule2;
    chemkit::Atom *O1 = molecule2.addAtom("O");
    chemkit::Atom *O2 = molecule2.addAtom("O");
    chemkit::Atom *O3 = molecule2.addAtom("O");
    chemkit::Bond *O1_O2 = molecule2.addBond(O1, O2);
    QVERIFY(molecule2.bond(O1, O2) == O1_O2);
    QVERIFY(molecule2.bond(C2, O1) == 0);
    QVERIFY(molecule2.bond(O3, C3) == 0);
    QVERIFY(molecule.bond(C1, O1) == 0);
    QVERIFY(molecule.bond(O1, C1) == 0);
    QVERIFY(molecule.bond(O3, C3) == 0);
}

void MoleculeTest::size()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.size(), size_t(0));

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), size_t(1));

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), size_t(2));

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), size_t(3));

    molecule.removeAtom(molecule.atoms()[0]);
    QCOMPARE(molecule.size(), size_t(2));

    molecule.clear();
    QCOMPARE(molecule.size(), size_t(0));
}

void MoleculeTest::isEmpty()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.isEmpty(), true);

    chemkit::Atom *atom = molecule.addAtom("C");
    QCOMPARE(molecule.isEmpty(), false);

    molecule.removeAtom(atom);
    QCOMPARE(molecule.isEmpty(), true);
}

void MoleculeTest::rings()
{
    chemkit::Molecule empty;
    QCOMPARE(size_t(empty.rings().size()), size_t(0));
    QCOMPARE(empty.ringCount(), size_t(0));

    chemkit::Molecule cyclopropane;
    chemkit::Atom *cyclopropane_C1 = cyclopropane.addAtom("C");
    chemkit::Atom *cyclopropane_C2 = cyclopropane.addAtom("C");
    chemkit::Atom *cyclopropane_C3 = cyclopropane.addAtom("C");
    QCOMPARE(size_t(cyclopropane.rings().size()), size_t(0));
    QCOMPARE(cyclopropane.ringCount(), size_t(0));

    cyclopropane.addBond(cyclopropane_C1, cyclopropane_C2);
    cyclopropane.addBond(cyclopropane_C2, cyclopropane_C3);
    QCOMPARE(size_t(cyclopropane.rings().size()), size_t(0));
    QCOMPARE(cyclopropane.ringCount(), size_t(0));

    cyclopropane.addBond(cyclopropane_C1, cyclopropane_C3);
    QCOMPARE(size_t(cyclopropane.rings().size()), size_t(1));
    QCOMPARE(cyclopropane.ringCount(), size_t(1));

    cyclopropane.removeAtom(cyclopropane_C2);
    QCOMPARE(size_t(cyclopropane.rings().size()), size_t(0));
    QCOMPARE(cyclopropane.ringCount(), size_t(0));
}

void MoleculeTest::distance()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    QCOMPARE(molecule.distance(H1, H2), chemkit::Real(0.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Real(0.0));

    H1->setPosition(1, 0, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Real(1.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Real(1.0));

    H2->setPosition(3, 0, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Real(2.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Real(2.0));

    H1->setPosition(0, 4, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Real(5.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Real(5.0));

    chemkit::Atom *H3 = molecule.addAtom("H");
    H3->setPosition(0, 0, -3);
    QCOMPARE(molecule.distance(H1, H3), chemkit::Real(5.0));
}

void MoleculeTest::center()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.center(), chemkit::Point3(0, 0, 0));

    chemkit::Atom *H1 = molecule.addAtom("H");
    QCOMPARE(molecule.center(), chemkit::Point3(0, 0, 0));

    H1->setPosition(1.0, 0, 0);
    QCOMPARE(molecule.center(), chemkit::Point3(1.0, 0, 0));

    chemkit::Atom *H2 = molecule.addAtom("H");
    H2->setPosition(3.0, 0, 0);
    QCOMPARE(molecule.center(), chemkit::Point3(2.0, 0, 0));
}

void MoleculeTest::bondAngle()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    C1->setPosition(1, 0, 0);
    chemkit::Atom *C2 = molecule.addAtom("C");
    C2->setPosition(0, 0, 0);
    chemkit::Atom *C3 = molecule.addAtom("C");
    C3->setPosition(0, 1, 0);
    QCOMPARE(qRound(molecule.bondAngle(C1, C2, C3)), 90);
    QCOMPARE(qRound(molecule.bondAngle(C3, C2, C1)), 90);

    C2->setPosition(0.5, 0.5, 0);
    QCOMPARE(qRound(molecule.bondAngle(C1, C2, C3)), 180);
    QCOMPARE(qRound(molecule.bondAngle(C3, C2, C1)), 180);

    C2->setPosition(0.25, 0.25, 0);
    QCOMPARE(qRound(molecule.bondAngle(C1, C2, C3)), 127);
    QCOMPARE(qRound(molecule.bondAngle(C3, C2, C1)), 127);
}

void MoleculeTest::torsionAngle()
{
}

void MoleculeTest::wilsonAngle()
{
}

void MoleculeTest::fragments()
{
    chemkit::Molecule neon;
    QCOMPARE(neon.fragmentCount(), size_t(0));

    neon.addAtom("Ne");
    QCOMPARE(neon.fragmentCount(), size_t(1));

    neon.addAtom("Ne");
    QCOMPARE(neon.fragmentCount(), size_t(2));

    neon.removeAtom(neon.atom(1));
    QCOMPARE(neon.fragmentCount(), size_t(1));

    neon.removeAtom(neon.atom(0));
    QCOMPARE(neon.fragmentCount(), size_t(0));
}

void MoleculeTest::isFragmented()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.isFragmented(), false);

    chemkit::Atom *C1 = molecule.addAtom("C");
    QCOMPARE(molecule.isFragmented(), false);

    chemkit::Atom *C2 = molecule.addAtom("C");
    QCOMPARE(molecule.isFragmented(), true);

    molecule.addBond(C1, C2);
    QCOMPARE(molecule.isFragmented(), false);

    chemkit::Atom *C3 = molecule.addAtom("C");
    QCOMPARE(molecule.isFragmented(), true);

    molecule.addBond(C2, C3);
    QCOMPARE(molecule.isFragmented(), false);

    molecule.removeBond(C1, C2);
    QCOMPARE(molecule.isFragmented(), true);

    molecule.clear();
    QCOMPARE(molecule.isFragmented(), false);
}

void MoleculeTest::removeFragment()
{
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    molecule.addBond(O1, H1);
    molecule.addBond(O1, H2);
    chemkit::Atom *O3 = molecule.addAtom("O");
    chemkit::Atom *H4 = molecule.addAtom("H");
    chemkit::Atom *H5 = molecule.addAtom("H");
    molecule.addBond(O3, H4);
    molecule.addBond(O3, H5);
    QCOMPARE(molecule.formula(), std::string("H4O2"));
    QCOMPARE(molecule.fragmentCount(), size_t(2));

    molecule.removeFragment(O3->fragment());
    QCOMPARE(molecule.formula(), std::string("H2O"));
    QCOMPARE(molecule.fragmentCount(), size_t(1));

    molecule.removeFragment(H1->fragment());
    QCOMPARE(molecule.formula(), std::string());
    QCOMPARE(molecule.fragmentCount(), size_t(0));
}

void MoleculeTest::rotate()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");

    C1->setPosition(0, -1, 0);
    C2->setPosition(0, 0, 0);
    C3->setPosition(1, 0, 0);

    molecule.coordinates()->rotate(chemkit::Vector3::UnitZ(), 90);
    QVERIFY(C1->position().isApprox(chemkit::Vector3(1, 0, 0)));
    QVERIFY(C2->position().isApprox(chemkit::Vector3(0, 0, 0)));
    QVERIFY(C3->position().isApprox(chemkit::Vector3(0, 1, 0)));
}

QTEST_APPLESS_MAIN(MoleculeTest)
