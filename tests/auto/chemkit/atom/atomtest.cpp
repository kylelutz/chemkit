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

#include "atomtest.h"

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void AtomTest::atomicNumber()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->atomicNumber(), 1);

    atom->setAtomicNumber(6);
    QCOMPARE(atom->atomicNumber(), 6);

    atom->setAtomicNumber(6);
    QCOMPARE(atom->atomicNumber(), 6);

    atom->setAtomicNumber(0);
    QCOMPARE(atom->atomicNumber(), 6);

    atom->setAtomicNumber(500);
    QCOMPARE(atom->atomicNumber(), 6);
}

void AtomTest::formalCharge()
{
    chemkit::Molecule molecule;
    chemkit::Atom *carbon = molecule.addAtom("C");
    QCOMPARE(carbon->formalCharge(), -4);

    for(int i = 0; i < 4; i++)
        molecule.addBond(carbon, molecule.addAtom("H"));
    QCOMPARE(carbon->formalCharge(), 0);

    for(int i = 0; i < 4; i++)
        molecule.addBond(carbon, molecule.addAtom("H"));
    QCOMPARE(carbon->formalCharge(), 4);
}

void AtomTest::partialCharge()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->partialCharge(), chemkit::Float(0.0));

    atom->setPartialCharge(2);
    QCOMPARE(atom->partialCharge(), chemkit::Float(2.0));
}

void AtomTest::symbol()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->symbol(), std::string("H"));

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(atom->symbol(), std::string("C"));
}

void AtomTest::name()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->name(), std::string("Hydrogen"));

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(atom->name(), std::string("Carbon"));
}

void AtomTest::electronegativity()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(qRound(atom->electronegativity()), 2);

    atom->setAtomicNumber(chemkit::Atom::Fluorine);
    QCOMPARE(qRound(atom->electronegativity()), 4);
}

void AtomTest::covalentRadius()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(qRound(atom->covalentRadius()), 0);

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(qRound(atom->covalentRadius()), 1);
}

void AtomTest::vanDerWaalsRadius()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(qRound(atom->vanDerWaalsRadius()), 1);

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(qRound(atom->vanDerWaalsRadius()), 2);
}

void AtomTest::is()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->is(chemkit::Atom::Hydrogen), true);
    QCOMPARE(atom->is(chemkit::Atom::Carbon), false);
    QCOMPARE(atom->isHeteroatom(), false);

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(atom->is(chemkit::Atom::Carbon), true);
    QCOMPARE(atom->is(chemkit::Atom::Oxygen), false);
    QCOMPARE(atom->isHeteroatom(), false);

    atom->setAtomicNumber(chemkit::Atom::Lithium);
    QCOMPARE(atom->is(chemkit::Atom::Lithium), true);
    QCOMPARE(atom->is(chemkit::Atom::Hydrogen), false);
    QCOMPARE(atom->isHeteroatom(), true);
}

void AtomTest::molecule()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QVERIFY(atom->molecule() == &molecule);
}

void AtomTest::rings()
{
    chemkit::Molecule benzene("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
    QCOMPARE(benzene.formula(), std::string("C6H6"));
    QCOMPARE(benzene.ringCount(), 1);
    chemkit::Ring *benzeneRing = benzene.rings()[0];

    foreach(chemkit::Atom *atom, benzene.atoms()){
        if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(atom->ringCount(), 0);
            QCOMPARE(atom->isInRing(), false);
            QCOMPARE(atom->isInRing(6), false);
            QCOMPARE(atom->isInRing(5), false);
            QCOMPARE(atom->isAromatic(), false);
            QVERIFY(atom->smallestRing() == 0);
        }
        else{
            QCOMPARE(atom->ringCount(), 1);
            QCOMPARE(atom->isInRing(), true);
            QCOMPARE(atom->isInRing(6), true);
            QCOMPARE(atom->isInRing(5), false);
            QCOMPARE(atom->isAromatic(), true);
            QVERIFY(atom->smallestRing() == benzeneRing);
        }
    }
}

void AtomTest::position()
{
    chemkit::Molecule molecule;
    chemkit::Atom *H1 = molecule.addAtom("H");
    QCOMPARE(H1->x(), chemkit::Float(0.0));
    QCOMPARE(H1->y(), chemkit::Float(0.0));
    QCOMPARE(H1->z(), chemkit::Float(0.0));

    H1->setPosition(4, 5, 6);
    QCOMPARE(H1->x(), chemkit::Float(4.0));
    QCOMPARE(H1->y(), chemkit::Float(5.0));
    QCOMPARE(H1->z(), chemkit::Float(6.0));

    H1->setPosition(chemkit::Point3());
    QCOMPARE(H1->x(), chemkit::Float(0.0));
    QCOMPARE(H1->y(), chemkit::Float(0.0));
    QCOMPARE(H1->z(), chemkit::Float(0.0));

    H1->setPosition(chemkit::Point3(-1, -2, -3));
    QCOMPARE(H1->x(), chemkit::Float(-1.0));
    QCOMPARE(H1->y(), chemkit::Float(-2.0));
    QCOMPARE(H1->z(), chemkit::Float(-3.0));
}

void AtomTest::distance()
{
    chemkit::Molecule molecule;
    chemkit::Atom *He1 = molecule.addAtom("He");
    chemkit::Atom *He2 = molecule.addAtom("He");
    QCOMPARE(He1->distance(He2), chemkit::Float(0.0));
    QCOMPARE(He2->distance(He1), chemkit::Float(0.0));
    QCOMPARE(He1->distance(He1), chemkit::Float(0.0));
    QCOMPARE(He2->distance(He2), chemkit::Float(0.0));

    He1->moveTo(0, 0, 5);
    QCOMPARE(He1->distance(He2), chemkit::Float(5.0));
    QCOMPARE(He2->distance(He1), chemkit::Float(5.0));
}

void AtomTest::pathTo()
{
    QList<chemkit::Atom *> atomPath;
    QList<chemkit::Bond *> bondPath;

    // propane
    chemkit::Molecule propane;
    chemkit::Atom *C1 = propane.addAtom("C");
    chemkit::Atom *C2 = propane.addAtom("C");
    chemkit::Atom *C3 = propane.addAtom("C");
    chemkit::Bond *C1_C2 = propane.addBond(C1, C2);
    chemkit::Bond *C2_C3 = propane.addBond(C2, C3);
    QCOMPARE(C1->atomCountTo(C2), 1);
    QCOMPARE(C1->atomCountTo(C3), 2);
    QCOMPARE(C1->atomCountTo(C1), 0);
    QCOMPARE(C1->bondCountTo(C2), 1);
    QCOMPARE(C1->bondCountTo(C3), 2);
    QCOMPARE(C3->bondCountTo(C1), 2);
    QCOMPARE(C1->bondCountTo(C1), 0);

    atomPath = C1->atomPathTo(C2);
    QCOMPARE(atomPath.size(), 1);
    QVERIFY(atomPath[0] == C2);

    atomPath = C1->atomPathTo(C3);
    QCOMPARE(atomPath.size(), 2);
    QVERIFY(atomPath[0] == C2);
    QVERIFY(atomPath[1] == C3);

    bondPath = C1->bondPathTo(C2);
    QCOMPARE(bondPath.size(), 1);
    QVERIFY(bondPath[0] == C1_C2);

    bondPath = C3->bondPathTo(C2);
    QCOMPARE(bondPath.size(), 1);
    QVERIFY(bondPath[0] == C2_C3);

    bondPath = C1->bondPathTo(C3);
    QCOMPARE(bondPath.size(), 2);
    QVERIFY(bondPath[0] == C1_C2);
    QVERIFY(bondPath[1] == C2_C3);

    QCOMPARE(C1->atomCountTo(C1, 0), 0);
    QCOMPARE(C1->atomCountTo(C1, 1), 0);
    QCOMPARE(C1->atomCountTo(C2, 0), 0);
    QCOMPARE(C1->atomCountTo(C2, 1), 1);
    QCOMPARE(C1->atomCountTo(C2, 2), 1);
    QCOMPARE(C1->atomCountTo(C3, 1), 0);
    QCOMPARE(C1->atomCountTo(C3, 2), 2);
    QCOMPARE(C1->atomCountTo(C3, 3), 2);

    QCOMPARE(C1->bondCountTo(C1, 0), 0);
    QCOMPARE(C1->bondCountTo(C1, 1), 0);
    QCOMPARE(C1->bondCountTo(C2, 0), 0);
    QCOMPARE(C1->bondCountTo(C2, 1), 1);
    QCOMPARE(C1->bondCountTo(C2, 2), 1);
    QCOMPARE(C1->bondCountTo(C3, 1), 0);
    QCOMPARE(C1->bondCountTo(C3, 2), 2);
    QCOMPARE(C1->bondCountTo(C3, 3), 2);

    // waters
    chemkit::Molecule waters;
    chemkit::Atom *O1 = waters.addAtom("O");
    chemkit::Atom *H2 = waters.addAtom("H");
    chemkit::Atom *H3 = waters.addAtom("H");
    chemkit::Atom *O4 = waters.addAtom("O");
    chemkit::Atom *H5 = waters.addAtom("H");
    chemkit::Atom *H6 = waters.addAtom("H");
    chemkit::Bond *O1_H2 = waters.addBond(O1, H2);
    chemkit::Bond *O1_H3 = waters.addBond(O1, H3);
    chemkit::Bond *O4_H5 = waters.addBond(O4, H5);
    chemkit::Bond *O4_H6 = waters.addBond(O4, H6);
    QCOMPARE(O1->atomCountTo(O4), 0);
    QCOMPARE(O1->bondCountTo(O4), 0);
    QCOMPARE(O1->atomCountTo(H2), 1);
    QCOMPARE(O4->atomCountTo(H6), 1);

    bondPath = H2->bondPathTo(H3);
    QCOMPARE(bondPath.size(), 2);
    QVERIFY(bondPath[0] == O1_H2);
    QVERIFY(bondPath[1] == O1_H3);

    bondPath = H6->bondPathTo(H5);
    QCOMPARE(bondPath.size(), 2);
    QVERIFY(bondPath[0] == O4_H6);
    QVERIFY(bondPath[1] == O4_H5);

    bondPath = H6->bondPathTo(H2);
    QCOMPARE(bondPath.size(), 0);

    // adenosine triphosphate
    chemkit::Molecule atp;
    O1 = atp.addAtom("O");
    C2 = atp.addAtom("C");
    C3 = atp.addAtom("C");
    chemkit::Atom *C4 = atp.addAtom("C");
    chemkit::Atom *C5 = atp.addAtom("C");
    chemkit::Atom *C6 = atp.addAtom("C");
    chemkit::Atom *O7 = atp.addAtom("O");
    chemkit::Atom *O8 = atp.addAtom("O");
    chemkit::Atom *N9 = atp.addAtom("N");
    chemkit::Atom *C10 = atp.addAtom("C");
    chemkit::Atom *N11 = atp.addAtom("N");
    chemkit::Atom *C12 = atp.addAtom("C");
    chemkit::Atom *C13 = atp.addAtom("C");
    chemkit::Atom *N14 = atp.addAtom("N");
    chemkit::Atom *C15 = atp.addAtom("C");
    chemkit::Atom *N16 = atp.addAtom("N");
    chemkit::Atom *C17 = atp.addAtom("C");
    chemkit::Atom *N18 = atp.addAtom("N");
    chemkit::Atom *O19 = atp.addAtom("O");
    chemkit::Atom *P20 = atp.addAtom("P");
    chemkit::Atom *O21 = atp.addAtom("O");
    chemkit::Atom *O22 = atp.addAtom("O");
    chemkit::Atom *O23 = atp.addAtom("O");
    chemkit::Atom *P24 = atp.addAtom("P");
    chemkit::Atom *O25 = atp.addAtom("O");
    chemkit::Atom *O26 = atp.addAtom("O");
    chemkit::Atom *O27 = atp.addAtom("O");
    chemkit::Atom *P28 = atp.addAtom("P");
    chemkit::Atom *O29 = atp.addAtom("O");
    chemkit::Atom *O30 = atp.addAtom("O");
    chemkit::Atom *O31 = atp.addAtom("O");
    chemkit::Bond *O1_C2 = atp.addBond(O1, C2);
    chemkit::Bond *O1_C5 = atp.addBond(O1, C5);
    C2_C3 = atp.addBond(C2, C3);
    chemkit::Bond *C2_N9 = atp.addBond(C2, N9);
    chemkit::Bond *C3_C4 = atp.addBond(C3, C4);
    chemkit::Bond *C3_O8 = atp.addBond(C3, O8);
    chemkit::Bond *C4_C5 = atp.addBond(C4, C5);
    chemkit::Bond *C4_O7 = atp.addBond(C4, O7);
    chemkit::Bond *C5_C6 = atp.addBond(C5, C6);
    chemkit::Bond *C6_O19 = atp.addBond(C6, O19);
    chemkit::Bond *N9_C10 = atp.addBond(N9, C10);
    chemkit::Bond *N9_C17 = atp.addBond(N9, C17);
    chemkit::Bond *C10_N11 = atp.addBond(C10, N11);
    chemkit::Bond *N11_C12 = atp.addBond(N11, C12);
    chemkit::Bond *C12_C13 = atp.addBond(C12, C13);
    chemkit::Bond *C12_C17 = atp.addBond(C12, C17);
    chemkit::Bond *C13_N14 = atp.addBond(C13, N14);
    chemkit::Bond *C13_N18 = atp.addBond(C13, N18);
    chemkit::Bond *N14_C15 = atp.addBond(N14, C15);
    chemkit::Bond *C15_N16 = atp.addBond(C15, N16);
    chemkit::Bond *N16_C17 = atp.addBond(N16, C17);
    chemkit::Bond *O19_P20 = atp.addBond(O19, P20);
    chemkit::Bond *P20_O21 = atp.addBond(P20, O21);
    chemkit::Bond *P20_O22 = atp.addBond(P20, O22);
    chemkit::Bond *P20_O23 = atp.addBond(P20, O23);
    chemkit::Bond *O23_P24 = atp.addBond(O23, P24);
    chemkit::Bond *P24_O25 = atp.addBond(P24, O25);
    chemkit::Bond *P24_O26 = atp.addBond(P24, O26);
    chemkit::Bond *P24_O27 = atp.addBond(P24, O27);
    chemkit::Bond *O27_P28 = atp.addBond(O27, P28);
    chemkit::Bond *P28_O29 = atp.addBond(P28, O29);
    chemkit::Bond *P28_O30 = atp.addBond(P28, O30);
    chemkit::Bond *P28_O31 = atp.addBond(P28, O31);
    QCOMPARE(O1->atomCountTo(O7), 3);
    QCOMPARE(O1->atomCountTo(O29), 9);
    QCOMPARE(O31->bondCountTo(N18), 15);

    atomPath = O1->atomPathTo(O23);
    QCOMPARE(atomPath.size(), 5);
    QVERIFY(atomPath[0] == C5);
    QVERIFY(atomPath[1] == C6);
    QVERIFY(atomPath[2] == O19);
    QVERIFY(atomPath[3] == P20);
    QVERIFY(atomPath[4] == O23);

    atomPath = N9->atomPathTo(N18);
    QCOMPARE(atomPath.size(), 4);
    QVERIFY(atomPath[0] == C17);
    QVERIFY(atomPath[1] == C12);
    QVERIFY(atomPath[2] == C13);
    QVERIFY(atomPath[3] == N18);

    bondPath = O7->bondPathTo(C15);
    QCOMPARE(bondPath.size(), 7);
    QVERIFY(bondPath[0] == C4_O7);
    QVERIFY(bondPath[1] == C3_C4);
    QVERIFY(bondPath[2] == C2_C3);
    QVERIFY(bondPath[3] == C2_N9);
    QVERIFY(bondPath[4] == N9_C17);
    QVERIFY(bondPath[5] == N16_C17);
    QVERIFY(bondPath[6] == C15_N16);

    bondPath = O21->bondPathTo(C13);
    QCOMPARE(bondPath.size(), 10);
    QVERIFY(bondPath[0] == P20_O21);
    QVERIFY(bondPath[1] == O19_P20);
    QVERIFY(bondPath[2] == C6_O19);
    QVERIFY(bondPath[3] == C5_C6);
    QVERIFY(bondPath[4] == O1_C5);
    QVERIFY(bondPath[5] == O1_C2);
    QVERIFY(bondPath[6] == C2_N9);
    QVERIFY(bondPath[7] == N9_C17);
    QVERIFY(bondPath[8] == C12_C17);
    QVERIFY(bondPath[9] == C12_C13);

    QCOMPARE(O29->atomCountTo(P24, 2), 0);
    QCOMPARE(O29->atomCountTo(P24, 3), 3);
    QCOMPARE(O29->atomCountTo(P24, 4), 3);
    QCOMPARE(O29->atomCountTo(O25, 3), 0);

    QCOMPARE(O29->bondCountTo(P24, 2), 0);
    QCOMPARE(O29->bondCountTo(P24, 3), 3);
    QCOMPARE(O29->bondCountTo(P24, 4), 3);
    QCOMPARE(O29->bondCountTo(O25, 3), 0);

    Q_UNUSED(C3_O8);
    Q_UNUSED(C4_C5);
    Q_UNUSED(N9_C10);
    Q_UNUSED(C10_N11);
    Q_UNUSED(N11_C12);
    Q_UNUSED(C13_N14);
    Q_UNUSED(C13_N18);
    Q_UNUSED(N14_C15);
    Q_UNUSED(P20_O22);
    Q_UNUSED(P20_O23);
    Q_UNUSED(O23_P24);
    Q_UNUSED(P24_O25);
    Q_UNUSED(P24_O26);
    Q_UNUSED(P24_O27);
    Q_UNUSED(O27_P28);
    Q_UNUSED(P28_O29);
    Q_UNUSED(P28_O30);
    Q_UNUSED(P28_O31);
}

void AtomTest::otherNeighbor()
{
    chemkit::Molecule water;
    chemkit::Atom *O1 = water.addAtom("O");
    chemkit::Atom *H2 = water.addAtom("H");
    chemkit::Atom *H3 = water.addAtom("H");
    water.addBond(O1, H2);
    water.addBond(O1, H3);
    QVERIFY(O1->otherNeighbor(H2) == H3);
    QVERIFY(O1->otherNeighbor(H3) == H2);
}

QTEST_APPLESS_MAIN(AtomTest)
