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

#include "bondtest.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void BondTest::atoms()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Bond *bond = molecule.addBond(H1, H2);
    QVERIFY(bond->atom1() == H1);
    QVERIFY(bond->atom2() == H2);
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
    QCOMPARE(bond->order(), chemkit::Bond::BondOrderType(1));

    bond->setOrder(chemkit::Bond::Double);
    QCOMPARE(bond->order(), chemkit::Bond::BondOrderType(2));

    bond->setOrder(chemkit::Bond::Triple);
    QCOMPARE(bond->order(), chemkit::Bond::BondOrderType(3));

    bond->setOrder(chemkit::Bond::Single);
    QCOMPARE(bond->order(), chemkit::Bond::BondOrderType(1));
}

void BondTest::is()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(C1, C2);
    QVERIFY(bond->is(chemkit::Bond::Single) == true);
    QVERIFY(bond->is(chemkit::Bond::Double) == false);

    bond->setOrder(chemkit::Bond::Double);
    QVERIFY(bond->is(chemkit::Bond::Double) == true);
    QVERIFY(bond->is(chemkit::Bond::Single) == false);
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
    QCOMPARE(benzene.formula(), std::string("C6H6"));
    QCOMPARE(benzene.ringCount(), size_t(1));
    chemkit::Ring *benzeneRing = benzene.rings()[0];

    foreach(chemkit::Bond *bond, benzene.bonds()){
        if(bond->contains(chemkit::Atom::Hydrogen)){
            QCOMPARE(bond->ringCount(), size_t(0));
            QCOMPARE(bond->isInRing(), false);
            QCOMPARE(bond->isInRing(6), false);
            QCOMPARE(bond->isInRing(5), false);
            QCOMPARE(bond->isAromatic(), false);
            QVERIFY(bond->smallestRing() == 0);
        }
        else{
            QCOMPARE(bond->ringCount(), size_t(1));
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
    QCOMPARE(bond->length(), chemkit::Real(0.0));

    H1->setPosition(0, 1, 0);
    QCOMPARE(bond->length(), chemkit::Real(1.0));

    H2->setPosition(0, -3, 0);
    QCOMPARE(bond->length(), chemkit::Real(4.0));
}

void BondTest::stereochemistry()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(C1, C2, chemkit::Bond::Double);
    QCOMPARE(bond->stereochemistry(), chemkit::Stereochemistry::None);

    bond->setStereochemistry(chemkit::Stereochemistry::Cis);
    QCOMPARE(bond->stereochemistry(), chemkit::Stereochemistry::Cis);

    bond->setStereochemistry(chemkit::Stereochemistry::Trans);
    QCOMPARE(bond->stereochemistry(), chemkit::Stereochemistry::Trans);

    bond->setStereochemistry(chemkit::Stereochemistry::None);
    QCOMPARE(bond->stereochemistry(), chemkit::Stereochemistry::None);
}

QTEST_APPLESS_MAIN(BondTest)
