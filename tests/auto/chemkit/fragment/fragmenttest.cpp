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

#include "fragmenttest.h"

#include <algorithm>

#include <chemkit/fragment.h>
#include <chemkit/molecule.h>

void FragmentTest::basic()
{
    chemkit::Molecule waters;
    chemkit::Atom *O1 = waters.addAtom("O");
    chemkit::Atom *H2 = waters.addAtom("H");
    chemkit::Atom *H3 = waters.addAtom("H");
    chemkit::Atom *O4 = waters.addAtom("O");
    chemkit::Atom *H5 = waters.addAtom("H");
    chemkit::Atom *H6 = waters.addAtom("H");
    waters.addBond(O1, H2);
    waters.addBond(O1, H3);
    waters.addBond(O4, H5);
    waters.addBond(O4, H6);
    QCOMPARE(waters.fragmentCount(), size_t(2));
}

void FragmentTest::molecule()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("H");
    const chemkit::Fragment *fragment = atom->fragment();
    QVERIFY(fragment->molecule() == &molecule);
}

void FragmentTest::atoms()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    molecule.addBond(C1, C2);

    const std::vector<chemkit::Atom *> C1_atoms = C1->fragment()->atoms();
    QCOMPARE(C1->fragment()->atomCount(), size_t(2));
    QCOMPARE(C1_atoms.size(), size_t(2));
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C1) != C1_atoms.end());
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C2) != C1_atoms.end());
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C3) == C1_atoms.end());

    for(size_t i = 0; i < C1_atoms.size(); i++){
        QVERIFY(C1->fragment()->atom(i) == C1_atoms[i]);
    }

    const std::vector<chemkit::Atom *> C3_atoms = C3->fragment()->atoms();
    QCOMPARE(C3->fragment()->atomCount(), size_t(1));
    QCOMPARE(C3_atoms.size(), size_t(1));
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C1) == C3_atoms.end());
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C2) == C3_atoms.end());
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C3) != C3_atoms.end());

    for(size_t i = 0; i < C3_atoms.size(); i++){
        QVERIFY(C3->fragment()->atom(i) == C3_atoms[i]);
    }
}

void FragmentTest::contains()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    QCOMPARE(C1->fragment()->contains(C1), true);
    QCOMPARE(C2->fragment()->contains(C1), false);
    QCOMPARE(C3->fragment()->contains(C3), true);

    molecule.addBond(C1, C2);
    QCOMPARE(C1->fragment()->contains(C1), true);
    QCOMPARE(C1->fragment()->contains(C2), true);
    QCOMPARE(C1->fragment()->contains(C3), false);
    QCOMPARE(C3->fragment()->contains(C1), false);
    QCOMPARE(C3->fragment()->contains(C2), false);
    QCOMPARE(C3->fragment()->contains(C3), true);
}

void FragmentTest::bonds()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    QCOMPARE(C1->fragment()->bondCount(), size_t(1));

    const std::vector<chemkit::Bond *> &C1_bonds = C1->fragment()->bonds();
    QVERIFY(std::find(C1_bonds.begin(), C1_bonds.end(), C1_C2) != C1_bonds.end());

    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    QCOMPARE(C2->fragment()->bondCount(), size_t(2));

    const std::vector<chemkit::Bond *> &C2_bonds = C2->fragment()->bonds();
    QVERIFY(std::find(C2_bonds.begin(), C2_bonds.end(), C2_C3) != C2_bonds.end());
}

QTEST_APPLESS_MAIN(FragmentTest)
