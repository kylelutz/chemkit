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

#include <QtTest>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/atommapping.h>

class MoleculeTest : public QObject
{
    Q_OBJECT

    private slots:
        void name();
        void formula();
        void mass();
        void addAtom();
        void addAtomCopy();
        void addBond();
        void bond();
        void size();
        void isEmpty();
        void substructure();
        void mapping();
        void find();
        void rings();
        void distance();
        void center();
        void bondAngle();
        void torsionAngle();
        void wilsonAngle();
        void fragments();
        void isFragmented();
        void removeFragment();
        void addConformer();
        void conformers();
};

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
    QCOMPARE(molecule.formula(), QString());

    molecule.addAtom("O");
    molecule.addAtom("O");
    QCOMPARE(molecule.formula(), QString("O2"));

    molecule.addAtom("Ac");
    QCOMPARE(molecule.formula(), QString("AcO2"));

    molecule.addAtom("C");
    QCOMPARE(molecule.formula(), QString("CAcO2"));

    molecule.addAtom("H");
    molecule.addAtom("H");
    QCOMPARE(molecule.formula(), QString("CH2AcO2"));

    molecule.clear();
    QCOMPARE(molecule.formula(), QString());
}

void MoleculeTest::mass()
{
    chemkit::Molecule molecule;
    QCOMPARE(molecule.mass(), chemkit::Float(0.0));

    molecule.addAtom("C");
    QCOMPARE(qRound(molecule.mass()), 12);

    molecule.addAtom("H");
    QCOMPARE(qRound(molecule.mass()), 13);

    molecule.clear();
    QCOMPARE(molecule.mass(), chemkit::Float(0.0));
}

void MoleculeTest::addAtom()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom1;
    chemkit::Atom *atom2;

    atom1 = molecule.addAtom(chemkit::Atom::Carbon);
    QVERIFY(atom1 != 0);
    QCOMPARE(atom1->atomicNumber(), 6);

    atom2 = molecule.addAtom(0);
    QVERIFY(atom2 == 0);

    atom2 = molecule.addAtom(-1);
    QVERIFY(atom2 == 0);

    atom2 = molecule.addAtom(500);
    QVERIFY(atom2 == 0);

    atom1 = molecule.addAtom("C");
    QVERIFY(atom1 != 0);
    QCOMPARE(atom1->symbol(), QString("C"));

    atom2 = molecule.addAtom(QString());
    QVERIFY(atom2 == 0);

    atom2 = molecule.addAtom("");
    QVERIFY(atom2 == 0);

    atom2 = molecule.addAtom("X");
    QVERIFY(atom2 == 0);
}

void MoleculeTest::addAtomCopy()
{
}

void MoleculeTest::addBond()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    QCOMPARE(molecule.bondCount(), 0);

    chemkit::Bond *H1_H2 = molecule.addBond(H1, H2);
    QVERIFY(H1_H2->contains(H1));
    QVERIFY(H1_H2->contains(H2));
    QCOMPARE(H1_H2->order(), 1);
    QCOMPARE(molecule.bondCount(), 1);
    QVERIFY(molecule.bonds()[0] == H1_H2);

    chemkit::Bond *H1_H2_again = molecule.addBond(H1, H2);
    QVERIFY(H1_H2_again == H1_H2);
    QCOMPARE(molecule.bondCount(), 1);

    chemkit::Bond *H1_H1 = molecule.addBond(H1, H1);
    QVERIFY(H1_H1 == 0);
    QCOMPARE(molecule.bondCount(), 1);
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
    QCOMPARE(molecule.size(), 0);

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), 1);

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), 2);

    molecule.addAtom("C");
    QCOMPARE(molecule.size(), 3);

    molecule.removeAtom(molecule.atoms()[0]);
    QCOMPARE(molecule.size(), 2);

    molecule.clear();
    QCOMPARE(molecule.size(), 0);
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

void MoleculeTest::substructure()
{
    chemkit::Molecule empty1;
    chemkit::Molecule empty2;
    QCOMPARE(empty1.isSubstructureOf(&empty2), true);
    QCOMPARE(empty1.contains(&empty2), true);
    QCOMPARE(empty2.isSubstructureOf(&empty1), true);
    QCOMPARE(empty2.contains(&empty1), true);

    chemkit::Molecule methane;
    methane.addAtom("C");

    chemkit::Molecule ethane;
    chemkit::Atom *ethane_C1 = ethane.addAtom("C");
    chemkit::Atom *ethane_C2 = ethane.addAtom("C");
    ethane.addBond(ethane_C1, ethane_C2);

    chemkit::Molecule propane;
    chemkit::Atom *propane_C1 = propane.addAtom("C");
    chemkit::Atom *propane_C2 = propane.addAtom("C");
    chemkit::Atom *propane_C3 = propane.addAtom("C");
    propane.addBond(propane_C1, propane_C2);
    propane.addBond(propane_C2, propane_C3);

    QCOMPARE(methane.isSubstructureOf(&methane), true);
    QCOMPARE(methane.isSubstructureOf(&ethane), true);
    QCOMPARE(methane.isSubstructureOf(&propane), true);
    QCOMPARE(methane.contains(&methane), true);
    QCOMPARE(methane.contains(&ethane), false);
    QCOMPARE(methane.contains(&propane), false);

    QCOMPARE(ethane.isSubstructureOf(&methane), false);
    QCOMPARE(ethane.isSubstructureOf(&ethane), true);
    QCOMPARE(ethane.isSubstructureOf(&propane), true);
    QCOMPARE(ethane.contains(&methane), true);
    QCOMPARE(ethane.contains(&ethane), true);
    QCOMPARE(ethane.contains(&propane), false);

    QCOMPARE(propane.isSubstructureOf(&methane), false);
    QCOMPARE(propane.isSubstructureOf(&ethane), false);
    QCOMPARE(propane.isSubstructureOf(&propane), true);
    QCOMPARE(propane.contains(&methane), true);
    QCOMPARE(propane.contains(&ethane), true);
    QCOMPARE(propane.contains(&propane), true);

    chemkit::Molecule benzene("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
    QCOMPARE(benzene.formula(), QString("C6H6"));
    QCOMPARE(benzene.isSubstructureOf(&benzene), true);

    QCOMPARE(methane.isSubstructureOf(&benzene), true);
    QCOMPARE(ethane.isSubstructureOf(&benzene), true);
    QCOMPARE(propane.isSubstructureOf(&benzene), false);

    chemkit::Molecule phenol("InChI=1/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H", "inchi");
    QCOMPARE(phenol.formula(), QString("C6H6O"));
    QCOMPARE(phenol.isSubstructureOf(&phenol), true);
    QCOMPARE(benzene.isSubstructureOf(&phenol), true);
    QCOMPARE(benzene.contains(&phenol), false);
    QCOMPARE(phenol.isSubstructureOf(&benzene), false);
    QCOMPARE(phenol.contains(&benzene), true);
}

void MoleculeTest::mapping()
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

    chemkit::AtomMapping mapping = methanol.mapping(&ethanol);
    QCOMPARE(mapping.size(), 2);
    QVERIFY(mapping.source() == &methanol);
    QVERIFY(mapping.target() == &ethanol);
    QVERIFY(mapping.map(ethanol_C2) == methanol_C1);
    QVERIFY(mapping.map(ethanol_O3) == methanol_O2);

    mapping = methanol.mapping(&ethanol, chemkit::Molecule::CompareHydrogens);
    QCOMPARE(mapping.size(), 3);
    QVERIFY(mapping.source() == &methanol);
    QVERIFY(mapping.target() == &ethanol);
    QVERIFY(mapping.map(ethanol_C2) == methanol_C1);
    QVERIFY(mapping.map(ethanol_O3) == methanol_O2);
    QVERIFY(mapping.map(ethanol_H4) == methanol_H3);
}

void MoleculeTest::find()
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

    chemkit::Molecule carboxyl;
    chemkit::Atom *carboxyl_C1 = carboxyl.addAtom("C");
    chemkit::Atom *carboxyl_O2 = carboxyl.addAtom("O");
    chemkit::Atom *carboxyl_O3 = carboxyl.addAtom("O");
    carboxyl.addBond(carboxyl_C1, carboxyl_O2, chemkit::Bond::Double);
    carboxyl.addBond(carboxyl_C1, carboxyl_O3);

    chemkit::Moiety carboxylMoiety = alanine.find(&carboxyl);
    QCOMPARE(carboxylMoiety.atomCount(), 3);
    QVERIFY(carboxylMoiety.atom(0) == alanine_C2);
    QVERIFY(carboxylMoiety.atom(1) == alanine_O5);
    QVERIFY(carboxylMoiety.atom(2) == alanine_O6);

    chemkit::Molecule carbonyl;
    chemkit::Atom *carbonyl_C1 = carbonyl.addAtom("C");
    chemkit::Atom *carbonyl_O2 = carbonyl.addAtom("O");
    carbonyl.addBond(carbonyl_C1, carbonyl_O2, chemkit::Bond::Double);

    chemkit::Moiety carbonylMoiety = carboxyl.find(&carbonyl);
    QCOMPARE(carbonylMoiety.atomCount(), 2);
    QVERIFY(carbonylMoiety.atom(0) == carboxyl_C1);
    QVERIFY(carbonylMoiety.atom(1) == carboxyl_O2);

    carboxylMoiety = carbonyl.find(&carboxyl);
    QCOMPARE(carboxylMoiety.isEmpty(), true);
}

void MoleculeTest::rings()
{
    chemkit::Molecule empty;
    QCOMPARE(empty.rings().size(), 0);
    QCOMPARE(empty.ringCount(), 0);

    chemkit::Molecule cyclopropane;
    chemkit::Atom *cyclopropane_C1 = cyclopropane.addAtom("C");
    chemkit::Atom *cyclopropane_C2 = cyclopropane.addAtom("C");
    chemkit::Atom *cyclopropane_C3 = cyclopropane.addAtom("C");
    QCOMPARE(cyclopropane.rings().size(), 0);
    QCOMPARE(cyclopropane.ringCount(), 0);

    cyclopropane.addBond(cyclopropane_C1, cyclopropane_C2);
    cyclopropane.addBond(cyclopropane_C2, cyclopropane_C3);
    QCOMPARE(cyclopropane.rings().size(), 0);
    QCOMPARE(cyclopropane.ringCount(), 0);

    cyclopropane.addBond(cyclopropane_C1, cyclopropane_C3);
    QCOMPARE(cyclopropane.rings().size(), 1);
    QCOMPARE(cyclopropane.ringCount(), 1);

    cyclopropane.removeAtom(cyclopropane_C2);
    QCOMPARE(cyclopropane.rings().size(), 0);
    QCOMPARE(cyclopropane.ringCount(), 0);
}

void MoleculeTest::distance()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    QCOMPARE(molecule.distance(H1, H2), chemkit::Float(0.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Float(0.0));

    H1->setPosition(1, 0, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Float(1.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Float(1.0));

    H2->setPosition(3, 0, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Float(2.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Float(2.0));

    H1->setPosition(0, 4, 0);
    QCOMPARE(molecule.distance(H1, H2), chemkit::Float(5.0));
    QCOMPARE(molecule.distance(H2, H1), chemkit::Float(5.0));

    chemkit::Atom *H3 = molecule.addAtom("H");
    H3->setPosition(0, 0, -3);
    QCOMPARE(molecule.distance(H1, H3), chemkit::Float(5.0));
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
    QCOMPARE(molecule.formula(), QString("H4O2"));
    QCOMPARE(molecule.fragmentCount(), 2);

    molecule.removeFragment(O3->fragment());
    QCOMPARE(molecule.formula(), QString("H2O"));
    QCOMPARE(molecule.fragmentCount(), 1);

    molecule.removeFragment(H1->fragment());
    QCOMPARE(molecule.formula(), QString());
    QCOMPARE(molecule.fragmentCount(), 0);
}

void MoleculeTest::addConformer()
{
    chemkit::Molecule molecule;
    chemkit::Atom *Na1 = molecule.addAtom("Na");
    chemkit::Atom *Cl2 = molecule.addAtom("Cl");
    QCOMPARE(molecule.center(), chemkit::Point3(0, 0, 0));

    chemkit::Conformer *conformer = molecule.addConformer();
    conformer->setPosition(Na1, chemkit::Point3(1, 0, 0));
    conformer->setPosition(Cl2, chemkit::Point3(3, 0, 0));
    QCOMPARE(molecule.center(), chemkit::Point3(0, 0, 0));

    chemkit::Conformer *originalConformer = molecule.conformer();
    molecule.setConformer(conformer);
    QCOMPARE(molecule.center(), chemkit::Point3(2, 0, 0));

    molecule.setConformer(originalConformer);
    QCOMPARE(molecule.center(), chemkit::Point3(0, 0, 0));
}

void MoleculeTest::conformers()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    QCOMPARE(molecule.conformerCount(), 1);
    QCOMPARE(molecule.conformers().size(), 1);
    QVERIFY(molecule.conformer() != 0);

    chemkit::Conformer *conformer = molecule.addConformer();
    QCOMPARE(molecule.conformerCount(), 2);
    conformer->setPosition(C1, chemkit::Point3(1, 2, 3));
    conformer->setPosition(C2, chemkit::Point3(2, 4, 6));

    QCOMPARE(C1->position(), chemkit::Point3(0, 0, 0));
    QCOMPARE(C2->position(), chemkit::Point3(0, 0, 0));

    molecule.setConformer(conformer);
    QCOMPARE(C1->position(), chemkit::Point3(1, 2, 3));
    QCOMPARE(C2->position(), chemkit::Point3(2, 4, 6));
}

QTEST_APPLESS_MAIN(MoleculeTest)
#include "moleculetest.moc"
