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

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void AtomTest::atomicNumber()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    QCOMPARE(atom->atomicNumber(), chemkit::Atom::AtomicNumberType(1));

    atom->setAtomicNumber(6);
    QCOMPARE(atom->atomicNumber(), chemkit::Atom::AtomicNumberType(6));

    atom->setAtomicNumber(6);
    QCOMPARE(atom->atomicNumber(), chemkit::Atom::AtomicNumberType(6));

    atom->setAtomicNumber(0);
    QCOMPARE(atom->atomicNumber(), chemkit::Atom::AtomicNumberType(6));

    atom->setAtomicNumber(200);
    QCOMPARE(atom->atomicNumber(), chemkit::Atom::AtomicNumberType(6));
}

void AtomTest::index()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    QCOMPARE(C1->index(), size_t(0));

    chemkit::Atom *C2 = molecule.addAtom("C");
    QCOMPARE(C2->index(), size_t(1));

    chemkit::Atom *C3 = molecule.addAtom("C");
    QCOMPARE(C3->index(), size_t(2));

    chemkit::Atom *C4 = molecule.addAtom("C");

    molecule.removeAtom(C2);
    QCOMPARE(C1->index(), size_t(0));
    QCOMPARE(C3->index(), size_t(1));
    QCOMPARE(C4->index(), size_t(2));

    molecule.removeAtom(C1);
    QCOMPARE(C3->index(), size_t(0));
    QCOMPARE(C4->index(), size_t(1));

    molecule.removeAtom(C3);
    QCOMPARE(C4->index(), size_t(0));
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
    QCOMPARE(atom->partialCharge(), chemkit::Real(0.0));

    atom->setPartialCharge(2);
    QCOMPARE(atom->partialCharge(), chemkit::Real(2.0));
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

    atom->setAtomicNumber(chemkit::Atom::Carbon);
    QCOMPARE(atom->is(chemkit::Atom::Carbon), true);
    QCOMPARE(atom->is(chemkit::Atom::Oxygen), false);

    atom->setAtomicNumber(chemkit::Atom::Lithium);
    QCOMPARE(atom->is(chemkit::Atom::Lithium), true);
    QCOMPARE(atom->is(chemkit::Atom::Hydrogen), false);
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
    QCOMPARE(benzene.ringCount(), size_t(1));
    chemkit::Ring *benzeneRing = benzene.rings()[0];

    foreach(chemkit::Atom *atom, benzene.atoms()){
        if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(atom->ringCount(), size_t(0));
            QCOMPARE(atom->isInRing(), false);
            QCOMPARE(atom->isInRing(6), false);
            QCOMPARE(atom->isInRing(5), false);
            QCOMPARE(atom->isAromatic(), false);
            QVERIFY(atom->smallestRing() == 0);
        }
        else{
            QCOMPARE(atom->ringCount(), size_t(1));
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
    QCOMPARE(H1->x(), chemkit::Real(0.0));
    QCOMPARE(H1->y(), chemkit::Real(0.0));
    QCOMPARE(H1->z(), chemkit::Real(0.0));

    H1->setPosition(4, 5, 6);
    QCOMPARE(H1->x(), chemkit::Real(4.0));
    QCOMPARE(H1->y(), chemkit::Real(5.0));
    QCOMPARE(H1->z(), chemkit::Real(6.0));

    H1->setPosition(chemkit::Point3(0, 0, 0));
    QCOMPARE(H1->x(), chemkit::Real(0.0));
    QCOMPARE(H1->y(), chemkit::Real(0.0));
    QCOMPARE(H1->z(), chemkit::Real(0.0));

    H1->setPosition(chemkit::Point3(-1, -2, -3));
    QCOMPARE(H1->x(), chemkit::Real(-1.0));
    QCOMPARE(H1->y(), chemkit::Real(-2.0));
    QCOMPARE(H1->z(), chemkit::Real(-3.0));
}

void AtomTest::distance()
{
    chemkit::Molecule molecule;
    chemkit::Atom *He1 = molecule.addAtom("He");
    chemkit::Atom *He2 = molecule.addAtom("He");
    QCOMPARE(He1->distance(He2), chemkit::Real(0.0));
    QCOMPARE(He2->distance(He1), chemkit::Real(0.0));
    QCOMPARE(He1->distance(He1), chemkit::Real(0.0));
    QCOMPARE(He2->distance(He2), chemkit::Real(0.0));

    He1->moveTo(0, 0, 5);
    QCOMPARE(He1->distance(He2), chemkit::Real(5.0));
    QCOMPARE(He2->distance(He1), chemkit::Real(5.0));
}

QTEST_APPLESS_MAIN(AtomTest)
