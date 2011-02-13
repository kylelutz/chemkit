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

class BondTest : public QObject
{
    Q_OBJECT

    private slots:
        void atoms();
        void otherAtom();
        void order();
        void molecule();
        void residue();
        void contains();
        void isTerminal();
        void rings();
        void polarity();
        void length();
};

void BondTest::atoms()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QVERIFY(bond->atom1() == H1);
    QVERIFY(bond->atom2() == H2);
    QVERIFY(bond->atoms().contains(H1));
    QVERIFY(bond->atoms().contains(H2));
}

void BondTest::otherAtom()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QVERIFY(bond->otherAtom(H1) == H2);
    QVERIFY(bond->otherAtom(H2) == H1);
}

void BondTest::order()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QCOMPARE(bond->order(), 1);

    bond->setOrder(chemkit::Bond::Double);
    QCOMPARE(bond->order(), 2);

    bond->setOrder(chemkit::Bond::Triple);
    QCOMPARE(bond->order(), 3);

    bond->setOrder(chemkit::Bond::Single);
    QCOMPARE(bond->order(), 1);
}

void BondTest::molecule()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QVERIFY(bond->molecule() == &molecule);
}

void BondTest::residue()
{
}

void BondTest::contains()
{
    chemkit::Molecule molecule;
    chemkit::Atom *hydrogen = molecule.addAtom("H");
    chemkit::Atom *carbon = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(hydrogen, carbon);
    QCOMPARE(bond->contains(hydrogen), true);
    QCOMPARE(bond->contains(carbon), true);
    QCOMPARE(bond->contains(chemkit::Atom::Hydrogen), true);
    QCOMPARE(bond->contains(chemkit::Atom::Carbon), true);
    QCOMPARE(bond->contains(chemkit::Atom::Oxygen), false);

    chemkit::Atom *oxygen = molecule.addAtom("O");
    QCOMPARE(bond->contains(oxygen), false);

    chemkit::Atom *nitrogen = molecule.addAtom("N");
    QCOMPARE(bond->contains(nitrogen), false);
    QCOMPARE(bond->containsBoth(hydrogen, carbon), true);
    QCOMPARE(bond->containsBoth(carbon, hydrogen), true);
    QCOMPARE(bond->containsBoth(carbon, carbon), true);
    QCOMPARE(bond->containsBoth(hydrogen, hydrogen), true);
    QCOMPARE(bond->containsBoth(chemkit::Atom::Hydrogen, chemkit::Atom::Carbon), true);
    QCOMPARE(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Hydrogen), true);
    QCOMPARE(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Carbon), false);
    QCOMPARE(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Oxygen), false);
    QCOMPARE(bond->containsBoth(carbon, oxygen), false);
    QCOMPARE(bond->containsBoth(carbon, nitrogen), false);
    QCOMPARE(bond->containsBoth(oxygen, nitrogen), false);
}

void BondTest::isTerminal()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    QCOMPARE(C1_C2->isTerminal(), true);
    QCOMPARE(C2_C3->isTerminal(), false);
    QCOMPARE(C3_C4->isTerminal(), true);
}

void BondTest::rings()
{
    chemkit::Molecule benzene("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
    QCOMPARE(benzene.formula(), QString("C6H6"));
    QCOMPARE(benzene.ringCount(), 1);
    chemkit::Ring *benzeneRing = benzene.rings()[0];

    foreach(chemkit::Bond *bond, benzene.bonds()){
        if(bond->contains(chemkit::Atom::Hydrogen)){
            QCOMPARE(bond->ringCount(), 0);
            QCOMPARE(bond->isInRing(), false);
            QCOMPARE(bond->isInRing(6), false);
            QCOMPARE(bond->isInRing(5), false);
            QCOMPARE(bond->isAromatic(), false);
            QVERIFY(bond->smallestRing() == 0);
        }
        else{
            QCOMPARE(bond->ringCount(), 1);
            QCOMPARE(bond->isInRing(), true);
            QCOMPARE(bond->isInRing(6), true);
            QCOMPARE(bond->isInRing(5), false);
            QCOMPARE(bond->isAromatic(), true);
            QVERIFY(bond->smallestRing() == benzeneRing);
        }
    }
}

void BondTest::polarity()
{
    chemkit::Molecule molecule;
    chemkit::Atom *hydrogen = molecule.addAtom("H");
    chemkit::Atom *carbon = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(hydrogen, carbon);
    QCOMPARE(qRound(bond->polarity()), 0);

    chemkit::Atom *oxygen = molecule.addAtom("O");
    bond = molecule.addBond(hydrogen, oxygen);
    QCOMPARE(qRound(bond->polarity()), 1);

    chemkit::Atom *sodium = molecule.addAtom("Na");
    chemkit::Atom *fluorine = molecule.addAtom("F");
    bond = molecule.addBond(sodium, fluorine);
    QCOMPARE(qRound(bond->polarity()), 3);

    chemkit::Atom *helium1 = molecule.addAtom("He");
    chemkit::Atom *helium2 = molecule.addAtom("He");
    bond = molecule.addBond(helium1, helium2);
    QCOMPARE(qRound(bond->polarity()), 0);
}

void BondTest::length()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QCOMPARE(bond->length(), chemkit::Float(0.0));

    H1->moveTo(0, 1, 0);
    QCOMPARE(bond->length(), chemkit::Float(1.0));

    H2->moveTo(0, -3, 0);
    QCOMPARE(bond->length(), chemkit::Float(4.0));
}

QTEST_APPLESS_MAIN(BondTest)
#include "bondtest.moc"
