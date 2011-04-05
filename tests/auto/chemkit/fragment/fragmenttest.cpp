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
    QCOMPARE(waters.fragmentCount(), 2);
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
    QCOMPARE(C1->fragment()->atomCount(), 2);
    QCOMPARE(C1_atoms.size(), 2UL);
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C1) != C1_atoms.end());
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C2) != C1_atoms.end());
    QVERIFY(std::find(C1_atoms.begin(), C1_atoms.end(), C3) == C1_atoms.end());

    const std::vector<chemkit::Atom *> C3_atoms = C3->fragment()->atoms();
    QCOMPARE(C3->fragment()->atomCount(), 1);
    QCOMPARE(C3_atoms.size(), 1UL);
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C1) == C3_atoms.end());
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C2) == C3_atoms.end());
    QVERIFY(std::find(C3_atoms.begin(), C3_atoms.end(), C3) != C3_atoms.end());
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
    QCOMPARE(C1->fragment()->bondCount(), 1);

    const std::vector<chemkit::Bond *> &C1_bonds = C1->fragment()->bonds();
    QVERIFY(std::find(C1_bonds.begin(), C1_bonds.end(), C1_C2) != C1_bonds.end());

    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    QCOMPARE(C2->fragment()->bondCount(), 2);

    const std::vector<chemkit::Bond *> &C2_bonds = C2->fragment()->bonds();
    QVERIFY(std::find(C2_bonds.begin(), C2_bonds.end(), C2_C3) != C2_bonds.end());
}

QTEST_APPLESS_MAIN(FragmentTest)
