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

// The ring-perception test verifies the ring perception and aromaticity
// perception algorithms against a large number of molecules. Each molecule
// is constructed and then each atom and bond it contains is checked for ring
// membership and aromaticity.

#include "ringperceptiontest.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>

void RingPerceptionTest::addHydrogens(chemkit::Molecule *molecule)
{
    std::vector<chemkit::Atom *> atoms(molecule->atoms().begin(),
                                       molecule->atoms().end());

    foreach(chemkit::Atom *atom, atoms){
        while(atom->valence() < atom->expectedValence()){
            molecule->addBond(atom, molecule->addAtom("H"));
        }
    }
}

/* anthracene (C14H10)
 *
 *      C1     C7    C11
 *    /   \\ /   \\ /   \\
 *  C6     C2     C8    C12
 *  ||      |      |     |
 *  C5     C3     C9    C13
 *    \   // \   // \   //
 *      C4     C10    C14
 */
void RingPerceptionTest::anthracene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *C14 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 2);
    chemkit::Bond *C3_C10 = molecule.addBond(C3, C10, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 1);
    chemkit::Bond *C8_C11 = molecule.addBond(C8, C11, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 2);
    chemkit::Bond *C9_C14 = molecule.addBond(C9, C14, 1);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 2);
    chemkit::Bond *C12_C13 = molecule.addBond(C12, C13, 1);
    chemkit::Bond *C13_C14 = molecule.addBond(C13, C14, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C14H10"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);
    chemkit::Ring *R3 = C11->smallestRing();
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);
    QCOMPARE(C13->isInRing(), true);
    QCOMPARE(C14->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);
    QCOMPARE(C13->isAromatic(), true);
    QCOMPARE(C14->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C8_C11->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C9_C14->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);
    QCOMPARE(C12_C13->isInRing(), true);
    QCOMPARE(C13_C14->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C10->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C8_C11->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C9_C14->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
    QCOMPARE(C12_C13->isAromatic(), true);
    QCOMPARE(C13_C14->isAromatic(), true);
}

/* anthraquinone (C14H8O2)
 *
 *            O15
 *            ||
 *      C1    C7     C11
 *   //   \  /   \  /   \\
 *  C6     C2     C8    C12
 *   |     ||     ||     |
 *  C5     C3     C9    C13
 *   \\   /  \   /  \   //
 *      C4    C10     C14
 *            ||
 *            O16
 */
void RingPerceptionTest::anthraquinone()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *C14 = molecule.addAtom("C");
    chemkit::Atom *O15 = molecule.addAtom("O");
    chemkit::Atom *O16 = molecule.addAtom("O");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_C10 = molecule.addBond(C3, C10, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 1);
    chemkit::Bond *C7_O15 = molecule.addBond(C7, O15, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 2);
    chemkit::Bond *C8_C11 = molecule.addBond(C8, C11, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    chemkit::Bond *C9_C14 = molecule.addBond(C9, C14, 1);
    chemkit::Bond *C10_O16 = molecule.addBond(C10, O16, 2);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 2);
    chemkit::Bond *C12_C13 = molecule.addBond(C12, C13, 1);
    chemkit::Bond *C13_C14 = molecule.addBond(C13, C14, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C14H8O2"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);
    chemkit::Ring *R3 = C11->smallestRing();
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);
    QCOMPARE(C13->isInRing(), true);
    QCOMPARE(C14->isInRing(), true);
    QCOMPARE(O15->isInRing(), false);
    QCOMPARE(O16->isInRing(), false);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);
    QCOMPARE(C13->isAromatic(), true);
    QCOMPARE(C14->isAromatic(), true);
    QCOMPARE(O15->isAromatic(), false);
    QCOMPARE(O16->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C7_O15->isInRing(), false);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C8_C11->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C9_C14->isInRing(), true);
    QCOMPARE(C10_O16->isInRing(), false);
    QCOMPARE(C11_C12->isInRing(), true);
    QCOMPARE(C12_C13->isInRing(), true);
    QCOMPARE(C13_C14->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C10->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C7_O15->isAromatic(), false);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C8_C11->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C9_C14->isAromatic(), true);
    QCOMPARE(C10_O16->isAromatic(), false);
    QCOMPARE(C11_C12->isAromatic(), true);
    QCOMPARE(C12_C13->isAromatic(), true);
    QCOMPARE(C13_C14->isAromatic(), true);
}

/* arsole (C4H5As)
 *
 *      As1
 *    /     \
 *  C5       C2
 *   \\     //
 *    C4 - C3
 */
void RingPerceptionTest::arsole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *As1 = molecule.addAtom("As");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *As1_C2 = molecule.addBond(As1, C2, 1);
    chemkit::Bond *As1_C5 = molecule.addBond(As1, C5, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H5As"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(As1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(As1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(As1_C2->isInRing(), true);
    QCOMPARE(As1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(As1_C2->isAromatic(), true);
    QCOMPARE(As1_C5->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* benzene (C6H6)
 *
 *      C1
 *   //    \
 *  C6     C2
 *   |     ||
 *  C5     C3
 *   \\    /
 *      C4
 */
void RingPerceptionTest::benzene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C6H6"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
}

/* benzimidazole (C7H6N2)
 *
 *      C1    N7
 *   //   \  /  \\
 *  C6     C2    C8
 *   |     ||    |
 *  C5     C3 -- N9
 *   \\   /
 *      C4
 */
void RingPerceptionTest::benzimidazole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *N7 = molecule.addAtom("N");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *N9 = molecule.addAtom("N");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C2_N7 = molecule.addBond(C2, N7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_N9 = molecule.addBond(C3, N9, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    chemkit::Bond *N7_C8 = molecule.addBond(N7, C8, 2);
    chemkit::Bond *C8_N9 = molecule.addBond(C8, N9, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C7H6N2"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = N7->smallestRing();
    QCOMPARE(R2->size(), size_t(5));
    QCOMPARE(R2->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(N7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(N9->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(N7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(N9->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_N7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_N9->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(N7_C8->isInRing(), true);
    QCOMPARE(C8_N9->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_N7->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_N9->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(N7_C8->isAromatic(), true);
    QCOMPARE(C8_N9->isAromatic(), true);
}

/* benzobicyclooctane (C12H14)
 *
 *      C2      C10
 *    /  |  \ //    \
 *  C5  C1   C9      C12
 *   |   |   |       ||
 *  C6  C3   C7      C11
 *    \  |  / \\    /
 *      C4       C8
 */
void RingPerceptionTest::benzobicyclooctane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C3 = molecule.addBond(C1, C3);
    chemkit::Bond *C2_C5 = molecule.addBond(C2, C5);
    chemkit::Bond *C2_C9 = molecule.addBond(C2, C9);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C6 = molecule.addBond(C4, C6);
    chemkit::Bond *C4_C7 = molecule.addBond(C4, C7);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C7_C9 = molecule.addBond(C7, C9);
    chemkit::Bond *C8_C11 = molecule.addBond(C8, C11);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 2);
    chemkit::Bond *C10_C12 = molecule.addBond(C10, C12);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C12H14"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R3 = C10->smallestRing();
    QVERIFY(R3 != 0);
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C3->isInRing(), true);
    QCOMPARE(C2_C5->isInRing(), true);
    QCOMPARE(C2_C9->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C6->isInRing(), true);
    QCOMPARE(C4_C7->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C7_C9->isInRing(), true);
    QCOMPARE(C8_C11->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C10_C12->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C3->isAromatic(), false);
    QCOMPARE(C2_C5->isAromatic(), false);
    QCOMPARE(C2_C9->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C6->isAromatic(), false);
    QCOMPARE(C4_C7->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C7_C9->isAromatic(), true);
    QCOMPARE(C8_C11->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C10_C12->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
}

/* benzonorborene (C11H12)
 *
 *      C3      C8
 *   /  |  \  /   \\
 *  C1  |   C5      C9
 *  |   C7  ||      |
 *  C2  |   C6      C10
 *   \  |  /  \   //
 *      C4      C11
 */
void RingPerceptionTest::benzonorborene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C3 = molecule.addBond(C1, C3, 1);
    chemkit::Bond *C2_C4 = molecule.addBond(C2, C4, 1);
    chemkit::Bond *C3_C5 = molecule.addBond(C3, C5, 1);
    chemkit::Bond *C3_C7 = molecule.addBond(C3, C7, 1);
    chemkit::Bond *C4_C6 = molecule.addBond(C4, C6, 1);
    chemkit::Bond *C4_C7 = molecule.addBond(C4, C7, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C5_C8 = molecule.addBond(C5, C8, 1);
    chemkit::Bond *C6_C11 = molecule.addBond(C6, C11, 1);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 2);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C11H12"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C5->smallestRing();
    QCOMPARE(R2->size(), size_t(5));
    QCOMPARE(R2->isAromatic(), false);
    chemkit::Ring *R3 = C8->smallestRing();
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C3->isInRing(), true);
    QCOMPARE(C2_C4->isInRing(), true);
    QCOMPARE(C3_C5->isInRing(), true);
    QCOMPARE(C3_C7->isInRing(), true);
    QCOMPARE(C4_C6->isInRing(), true);
    QCOMPARE(C4_C7->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C8->isInRing(), true);
    QCOMPARE(C6_C11->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C3->isAromatic(), false);
    QCOMPARE(C2_C4->isAromatic(), false);
    QCOMPARE(C3_C5->isAromatic(), false);
    QCOMPARE(C3_C7->isAromatic(), false);
    QCOMPARE(C4_C6->isAromatic(), false);
    QCOMPARE(C4_C7->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C5_C8->isAromatic(), true);
    QCOMPARE(C6_C11->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C10_C11->isAromatic(), true);
}

/* bicyclooctane (C8H14)
 *
 *      C2
 *    /  |  \
 *  C1  C7  C3
 *   |   |   |
 *  C6  C8  C4
 *    \  |  /
 *       C5
 */
void RingPerceptionTest::bicyclooctane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C5_C8 = molecule.addBond(C5, C8);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H14"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C3->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C8->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C6->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C2_C7->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C5_C8->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
}

/* biotin (C10H16N2O3S)
 *
 *       O1
 *       ||
 *       C2
 *    /      \
 *   N9      N3
 *    \      /
 *    C8 -- C4
 *    /      \
 *   C7       C5      C11     C13
 *    \      /   \   /   \   /   \
 *       S6       C10     C12    C14 -- O16
 *                               ||
 *                               O15
 */
void RingPerceptionTest::biotin()
{
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *N3 = molecule.addAtom("N");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *S6 = molecule.addAtom("S");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *N9 = molecule.addAtom("N");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *C14 = molecule.addAtom("C");
    chemkit::Atom *O15 = molecule.addAtom("O");
    chemkit::Atom *O16 = molecule.addAtom("O");
    chemkit::Bond *O1_C2 = molecule.addBond(O1, C2, 2);
    chemkit::Bond *C2_N3 = molecule.addBond(C2, N3, 1);
    chemkit::Bond *C2_N9 = molecule.addBond(C2, N9, 1);
    chemkit::Bond *N3_C4 = molecule.addBond(N3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C4_C8 = molecule.addBond(C4, C8, 1);
    chemkit::Bond *C5_S6 = molecule.addBond(C5, S6, 1);
    chemkit::Bond *C5_C10 = molecule.addBond(C5, C10, 1);
    chemkit::Bond *S6_C7 = molecule.addBond(S6, C7, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 1);
    chemkit::Bond *C8_N9 = molecule.addBond(C8, N9, 1);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 1);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 1);
    chemkit::Bond *C12_C13 = molecule.addBond(C12, C13, 1);
    chemkit::Bond *C13_C14 = molecule.addBond(C13, C14, 1);
    chemkit::Bond *C14_O15 = molecule.addBond(C14, O15, 2);
    chemkit::Bond *C14_O16 = molecule.addBond(C14, O16, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C10H16N2O3S"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C2->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = S6->smallestRing();
    QCOMPARE(R2->size(), size_t(5));
    QCOMPARE(R2->isAromatic(), false);

    QCOMPARE(O1->isInRing(), false);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(N3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(S6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(N9->isInRing(), true);
    QCOMPARE(C10->isInRing(), false);
    QCOMPARE(C11->isInRing(), false);
    QCOMPARE(C12->isInRing(), false);
    QCOMPARE(C13->isInRing(), false);
    QCOMPARE(C14->isInRing(), false);
    QCOMPARE(O15->isInRing(), false);
    QCOMPARE(O16->isInRing(), false);

    QCOMPARE(O1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(N3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(S6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(N9->isAromatic(), false);
    QCOMPARE(C10->isAromatic(), false);
    QCOMPARE(C11->isAromatic(), false);
    QCOMPARE(C12->isAromatic(), false);
    QCOMPARE(C13->isAromatic(), false);
    QCOMPARE(C14->isAromatic(), false);
    QCOMPARE(O15->isAromatic(), false);
    QCOMPARE(O16->isAromatic(), false);

    QCOMPARE(O1_C2->isInRing(), false);
    QCOMPARE(C2_N3->isInRing(), true);
    QCOMPARE(C2_N9->isInRing(), true);
    QCOMPARE(N3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C4_C8->isInRing(), true);
    QCOMPARE(C5_S6->isInRing(), true);
    QCOMPARE(C5_C10->isInRing(), false);
    QCOMPARE(S6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_N9->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), false);
    QCOMPARE(C11_C12->isInRing(), false);
    QCOMPARE(C12_C13->isInRing(), false);
    QCOMPARE(C13_C14->isInRing(), false);
    QCOMPARE(C14_O15->isInRing(), false);
    QCOMPARE(C14_O16->isInRing(), false);

    QCOMPARE(O1_C2->isAromatic(), false);
    QCOMPARE(C2_N3->isAromatic(), false);
    QCOMPARE(C2_N9->isAromatic(), false);
    QCOMPARE(N3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C4_C8->isAromatic(), false);
    QCOMPARE(C5_S6->isAromatic(), false);
    QCOMPARE(C5_C10->isAromatic(), false);
    QCOMPARE(S6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C8_N9->isAromatic(), false);
    QCOMPARE(C10_C11->isAromatic(), false);
    QCOMPARE(C11_C12->isAromatic(), false);
    QCOMPARE(C12_C13->isAromatic(), false);
    QCOMPARE(C13_C14->isAromatic(), false);
    QCOMPARE(C14_O15->isAromatic(), false);
    QCOMPARE(C14_O16->isAromatic(), false);
}

/* biphenylene (C12H8)
 *
 *      C1           C9
 *   //   \        /   \\
 *  C6     C2 -- C7     C10
 *   |     ||    ||      |
 *  C5     C3 -- C8     C11
 *   \\   /        \   //
 *      C4           C12
 */
void RingPerceptionTest::biphenylene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_C8 = molecule.addBond(C3, C8, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C7_C9 = molecule.addBond(C7, C9, 1);
    chemkit::Bond *C8_C12 = molecule.addBond(C8, C12, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 2);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 1);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C12H8"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C2->smallestRing();
    QCOMPARE(R2->size(), size_t(4));
    QCOMPARE(R2->isAromatic(), false);
    chemkit::Ring *R3 = C9->smallestRing();
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C8->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C7_C9->isInRing(), true);
    QCOMPARE(C8_C12->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C8->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C7_C9->isAromatic(), true);
    QCOMPARE(C8_C12->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C10_C11->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
}

/* cubane (C8H8)
 *
 *  C8  ------ C5
 *  | \       / |
 *  |  C1 - C4  |
 *  |  |     |  |
 *  |  C2 - C3  |
 *  | /       \ |
 *  C7 ------- C6
 */
void RingPerceptionTest::cubane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C4 = molecule.addBond(C1, C4);
    chemkit::Bond *C1_C8 = molecule.addBond(C1, C8);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C6 = molecule.addBond(C3, C6);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C5_C8 = molecule.addBond(C5, C8);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H8"));

    QCOMPARE(molecule.ringCount(), size_t(5));

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C4->isInRing(), true);
    QCOMPARE(C1_C8->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C6->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C8->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C4->isAromatic(), false);
    QCOMPARE(C1_C8->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C2_C7->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C6->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C5_C8->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
}

/* cyclobutane (C4H8)
 *
 *  C1 - C2
 *  |     |
 *  C4 - C3
 */
void RingPerceptionTest::cyclobutane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C4 = molecule.addBond(C1, C4);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H8"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(4));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C4->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C4->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
}

/* cyclodecane (C10H20)
 *
 *     C1        C3
 *   /    \    /    \
 *  C10     C2       C4
 *  |                 |
 *  C9      C7       C5
 *   \    /    \    /
 *     C8        C6
 */
void RingPerceptionTest::cyclodecane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C10 = molecule.addBond(C1, C10);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C10H20"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(10));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(C9->isAromatic(), false);
    QCOMPARE(C10->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C10->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C10->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C8_C9->isAromatic(), false);
    QCOMPARE(C9_C10->isAromatic(), false);
}

/* cycloheptane (C7H14)
 *
 *        C1
 *     /      \
 *    C7      C2
 *    |        |
 *    C6      C3
 *     \      /
 *     C5 -- C4
 */
void RingPerceptionTest::cycloheptane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C7 = molecule.addBond(C1, C7);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C7H14"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(7));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C7->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C7->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
}

/* cyclohexane (C6H12)
 *
 *      C1
 *    /    \
 *  C6     C2
 *   |      |
 *  C5     C3
 *    \    /
 *      C4
 */
void RingPerceptionTest::cyclohexane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C6H12"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C6->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
}

/* cyclononane (C9H18)
 *
 *      C1
 *   /      \
 *  C9      C2
 *  |        |
 *  C8      C3
 *  |        |
 *  C7      C4
 *   \      /
 *    C6 - C5
 */
void RingPerceptionTest::cyclononane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C9 = molecule.addBond(C1, C9);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C9H18"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(9));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(C9->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C9->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C9->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C8_C9->isAromatic(), false);
}

/* cyclooctane (C8H16)
 *
 *    C1 -- C2
 *   /        \
 *  C8        C3
 *  |          |
 *  C7        C4
 *   \        /
 *    C6 -- C5
 */
void RingPerceptionTest::cyclooctane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C8 = molecule.addBond(C1, C8);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H16"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(8));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C8->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C8->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
}

/* cyclooctatetraene (C8H8)
 *
 *    C1 == C2
 *   /        \
 *  C8        C3
 *  ||        ||
 *  C7        C4
 *   \        /
 *    C6 == C5
 */
void RingPerceptionTest::cyclooctatetraene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 2);
    chemkit::Bond *C1_C8 = molecule.addBond(C1, C8, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 2);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H8"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(8));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C8->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C8->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
}

/* cyclopentane (C5H10)
 *
 *     C1
 *   /    \
 *  C5    C2
 *  |      |
 *  C4 -- C3
 */
void RingPerceptionTest::cyclopentane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C5 = molecule.addBond(C1, C5);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C5H10"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C5->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
}

/* cyclopropane (C3H6)
 *
 *     C1
 *   /    \
 *  C3 -- C2
 */
void RingPerceptionTest::cyclopropane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C3 = molecule.addBond(C1, C3);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C3H6"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(3));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C3->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C3->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
}

/* decalin (C10H18)
 *
 *      C1     C7
 *    /    \ /    \
 *  C6     C2     C8
 *   |      |      |
 *  C5     C3     C9
 *    \    / \    /
 *      C4     C10
 */
void RingPerceptionTest::decalin()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C10 = molecule.addBond(C10, C3);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C10H18"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(C9->isAromatic(), false);
    QCOMPARE(C10->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C6->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C2_C7->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C10->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C8_C9->isAromatic(), false);
    QCOMPARE(C9_C10->isAromatic(), false);
}

/* furan (C4H4O)
 *
 *      O1
 *    /    \
 *  C5      C2
 *   \\    //
 *   C4 - C3
 */
void RingPerceptionTest::furan()
{
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *O1_C2 = molecule.addBond(O1, C2, 1);
    chemkit::Bond *O1_C5 = molecule.addBond(C5, O1, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H4O"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(O1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(O1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(O1_C2->isInRing(), true);
    QCOMPARE(O1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(O1_C2->isAromatic(), true);
    QCOMPARE(O1_C5->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* indole (C8H7N)
 *
 *     N1     C6
 *   /   \  /   \\
 *  C5    C2     C7
 *  ||    ||      |
 *  C4 -- C3     C8
 *          \   //
 *            C9
 */
void RingPerceptionTest::indole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Bond *N1_C2 = molecule.addBond(N1, C2, 1);
    chemkit::Bond *N1_C5 = molecule.addBond(N1, C5, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C2_C6 = molecule.addBond(C2, C6, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_C9 = molecule.addBond(C3, C9, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7, 2);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 1);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H7N"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = N1->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C6->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);

    QCOMPARE(N1_C2->isInRing(), true);
    QCOMPARE(N1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C6->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C9->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);

    QCOMPARE(N1_C2->isAromatic(), true);
    QCOMPARE(N1_C5->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C6->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C9->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C6_C7->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
}

/* imidazole (C3H4N2)
 *
 *      N1
 *    /    \
 *  C5      C2
 *   \\    //
 *   C4 - N3
 */
void RingPerceptionTest::imidazole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *N3 = molecule.addAtom("N");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *N1_C2 = molecule.addBond(N1, C2, 1);
    chemkit::Bond *N1_C5 = molecule.addBond(N1, C5, 1);
    chemkit::Bond *C2_N3 = molecule.addBond(C2, N3, 2);
    chemkit::Bond *N3_C4 = molecule.addBond(N3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C3H4N2"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(N3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(N3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(N1_C2->isInRing(), true);
    QCOMPARE(N1_C5->isInRing(), true);
    QCOMPARE(C2_N3->isInRing(), true);
    QCOMPARE(N3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(N1_C2->isAromatic(), true);
    QCOMPARE(N1_C5->isAromatic(), true);
    QCOMPARE(C2_N3->isAromatic(), true);
    QCOMPARE(N3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* ladderane (C12H16)
 *
 *  C1 - C7 - C5 - C3 - C9  - C11
 *  |    |    |    |    |     |
 *  C2 - C8 - C6 - C4 - C10 - C12
 */
void RingPerceptionTest::ladderane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C7 = molecule.addBond(C1, C7);
    chemkit::Bond *C2_C8 = molecule.addBond(C2, C8);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C5 = molecule.addBond(C3, C5);
    chemkit::Bond *C3_C9 = molecule.addBond(C3, C9);
    chemkit::Bond *C4_C6 = molecule.addBond(C4, C6);
    chemkit::Bond *C4_C10 = molecule.addBond(C4, C10);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C5_C7 = molecule.addBond(C5, C7);
    chemkit::Bond *C6_C8 = molecule.addBond(C6, C8);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10);
    chemkit::Bond *C9_C11 = molecule.addBond(C9, C11);
    chemkit::Bond *C10_C12 = molecule.addBond(C10, C12);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C12H16"));

    QCOMPARE(molecule.ringCount(), size_t(5));
    foreach(const chemkit::Ring *ring, molecule.rings())
        QCOMPARE(ring->size(), size_t(4));

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(C9->isAromatic(), false);
    QCOMPARE(C10->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C7->isInRing(), true);
    QCOMPARE(C2_C8->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C5->isInRing(), true);
    QCOMPARE(C3_C9->isInRing(), true);
    QCOMPARE(C4_C6->isInRing(), true);
    QCOMPARE(C4_C10->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C7->isInRing(), true);
    QCOMPARE(C6_C8->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C9_C11->isInRing(), true);
    QCOMPARE(C10_C12->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C7->isAromatic(), false);
    QCOMPARE(C2_C8->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C5->isAromatic(), false);
    QCOMPARE(C3_C9->isAromatic(), false);
    QCOMPARE(C4_C6->isAromatic(), false);
    QCOMPARE(C4_C10->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C5_C7->isAromatic(), false);
    QCOMPARE(C6_C8->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C9_C10->isAromatic(), false);
    QCOMPARE(C9_C11->isAromatic(), false);
    QCOMPARE(C10_C12->isAromatic(), false);
    QCOMPARE(C11_C12->isAromatic(), false);
}

/* naphthalene (C10H8)
 *
 *      C1     C7
 *    /   \\ /   \\
 *  C6     C2     C8
 *  ||      |      |
 *  C5     C3     C9
 *    \   // \   //
 *      C4     C10
 */
void RingPerceptionTest::naphthalene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 2);
    chemkit::Bond *C3_C10 = molecule.addBond(C3, C10, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C10H8"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C10->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
}

/* norbornane (C7H12)
 *
 *     C1 - C2
 *    /       \
 *   C6 - C7 - C3
 *    \       /
 *     C5 - C4
 */
void RingPerceptionTest::norbornane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C7 = molecule.addBond(C3, C7);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C7H12"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C4->smallestRing();
    QCOMPARE(R2->size(), size_t(5));
    QCOMPARE(R2->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C7->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C6->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C7->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
}

/* oxazole (C3H3NO)
 *
 *      O1
 *    /    \
 *  C5      C2
 *   \\    //
 *   C4 - N3
 */
void RingPerceptionTest::oxazole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *N3 = molecule.addAtom("N");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *O1_C2 = molecule.addBond(O1, C2, 1);
    chemkit::Bond *O1_C5 = molecule.addBond(O1, C5, 1);
    chemkit::Bond *C2_N3 = molecule.addBond(C2, N3, 2);
    chemkit::Bond *N3_C4 = molecule.addBond(N3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C3H3NO"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(O1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(N3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(O1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(N3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(O1_C2->isInRing(), true);
    QCOMPARE(O1_C5->isInRing(), true);
    QCOMPARE(C2_N3->isInRing(), true);
    QCOMPARE(N3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(O1_C2->isAromatic(), true);
    QCOMPARE(O1_C5->isAromatic(), true);
    QCOMPARE(C2_N3->isAromatic(), true);
    QCOMPARE(N3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* oxirane (C2H4O)
 *
 *     O1
 *   /    \
 *  C3 -- C2
 */
void RingPerceptionTest::oxirane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Bond *O1_C2 = molecule.addBond(O1, C2);
    chemkit::Bond *O1_C3 = molecule.addBond(O1, C3);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C2H4O"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(3));
    QCOMPARE(R1->isAromatic(), false);

    QCOMPARE(O1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);

    QCOMPARE(O1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);

    QCOMPARE(O1_C2->isInRing(), true);
    QCOMPARE(O1_C3->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);

    QCOMPARE(O1_C2->isAromatic(), false);
    QCOMPARE(O1_C3->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
}

/* porphin (C20H14N4)
 *
 *      C1      C3      C5
 *   //    \  //   \  /   \\
 *  C23     C2      C4     C6
 *   \     /        \\     /
 *  C22 - N24        N8 - C7
 *   //                   \\
 *  C21                    C9
 *   \                     /
 *  C19 = N20       N14 - C10
 *   /      \        /     \\
 *  C18    C16     C13     C11
 *   \\    /  \\   /  \\   /
 *     C17      C15     C12
 */
void RingPerceptionTest::porphin()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *N8 = molecule.addAtom("N");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *N14 = molecule.addAtom("N");
    chemkit::Atom *C15 = molecule.addAtom("C");
    chemkit::Atom *C16 = molecule.addAtom("C");
    chemkit::Atom *C17 = molecule.addAtom("C");
    chemkit::Atom *C18 = molecule.addAtom("C");
    chemkit::Atom *C19 = molecule.addAtom("C");
    chemkit::Atom *N20 = molecule.addAtom("N");
    chemkit::Atom *C21 = molecule.addAtom("C");
    chemkit::Atom *C22 = molecule.addAtom("C");
    chemkit::Atom *C23 = molecule.addAtom("C");
    chemkit::Atom *N24 = molecule.addAtom("N");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C23 = molecule.addBond(C1, C23, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C2_N24 = molecule.addBond(C2, N24, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C4_N8 = molecule.addBond(C4, N8, 2);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7, 1);
    chemkit::Bond *C7_N8 = molecule.addBond(C7, N8, 1);
    chemkit::Bond *C7_C9 = molecule.addBond(C7, C9, 2);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 2);
    chemkit::Bond *C10_N14 = molecule.addBond(C10, N14, 1);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 1);
    chemkit::Bond *C12_C13 = molecule.addBond(C12, C13, 2);
    chemkit::Bond *C13_N14 = molecule.addBond(C13, N14, 1);
    chemkit::Bond *C13_C15 = molecule.addBond(C13, C15, 1);
    chemkit::Bond *C15_C16 = molecule.addBond(C15, C16, 2);
    chemkit::Bond *C16_C17 = molecule.addBond(C16, C17, 1);
    chemkit::Bond *C16_N20 = molecule.addBond(C16, N20, 1);
    chemkit::Bond *C17_C18 = molecule.addBond(C17, C18, 2);
    chemkit::Bond *C18_C19 = molecule.addBond(C18, C19, 1);
    chemkit::Bond *C19_N20 = molecule.addBond(C19, N20, 2);
    chemkit::Bond *C19_C21 = molecule.addBond(C19, C21, 1);
    chemkit::Bond *C21_C22 = molecule.addBond(C21, C22, 2);
    chemkit::Bond *C22_C23 = molecule.addBond(C22, C23, 1);
    chemkit::Bond *C22_N24 = molecule.addBond(C22, N24, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C20H14N4"));

    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C3->smallestRing();
    QCOMPARE(R2->size(), size_t(16));
    QCOMPARE(R2->isAromatic(), true);
    chemkit::Ring *R3 = C4->smallestRing();
    QCOMPARE(R3->size(), size_t(5));
    QCOMPARE(R3->isAromatic(), true);
    chemkit::Ring *R4 = C10->smallestRing();
    QCOMPARE(R4->size(), size_t(5));
    QCOMPARE(R4->isAromatic(), true);
    chemkit::Ring *R5 = C16->smallestRing();
    QCOMPARE(R5->size(), size_t(5));
    QCOMPARE(R5->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(N8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);
    QCOMPARE(C13->isInRing(), true);
    QCOMPARE(N14->isInRing(), true);
    QCOMPARE(C15->isInRing(), true);
    QCOMPARE(C16->isInRing(), true);
    QCOMPARE(C17->isInRing(), true);
    QCOMPARE(C18->isInRing(), true);
    QCOMPARE(C19->isInRing(), true);
    QCOMPARE(N20->isInRing(), true);
    QCOMPARE(C21->isInRing(), true);
    QCOMPARE(C22->isInRing(), true);
    QCOMPARE(C23->isInRing(), true);
    QCOMPARE(N24->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(N8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);
    QCOMPARE(C13->isAromatic(), true);
    QCOMPARE(N14->isAromatic(), true);
    QCOMPARE(C15->isAromatic(), true);
    QCOMPARE(C16->isAromatic(), true);
    QCOMPARE(C17->isAromatic(), true);
    QCOMPARE(C18->isAromatic(), true);
    QCOMPARE(C19->isAromatic(), true);
    QCOMPARE(N20->isAromatic(), true);
    QCOMPARE(C21->isAromatic(), true);
    QCOMPARE(C22->isAromatic(), true);
    QCOMPARE(C23->isAromatic(), true);
    QCOMPARE(N24->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C23->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_N24->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C4_N8->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_N8->isInRing(), true);
    QCOMPARE(C7_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), true);
    QCOMPARE(C10_N14->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);
    QCOMPARE(C12_C13->isInRing(), true);
    QCOMPARE(C13_N14->isInRing(), true);
    QCOMPARE(C13_C15->isInRing(), true);
    QCOMPARE(C15_C16->isInRing(), true);
    QCOMPARE(C16_C17->isInRing(), true);
    QCOMPARE(C16_N20->isInRing(), true);
    QCOMPARE(C17_C18->isInRing(), true);
    QCOMPARE(C18_C19->isInRing(), true);
    QCOMPARE(C19_N20->isInRing(), true);
    QCOMPARE(C21_C22->isInRing(), true);
    QCOMPARE(C22_C23->isInRing(), true);
    QCOMPARE(C22_N24->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C23->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_N24->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C4_N8->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C6_C7->isAromatic(), true);
    QCOMPARE(C7_N8->isAromatic(), true);
    QCOMPARE(C7_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C10_C11->isAromatic(), true);
    QCOMPARE(C10_N14->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
    QCOMPARE(C12_C13->isAromatic(), true);
    QCOMPARE(C13_N14->isAromatic(), true);
    QCOMPARE(C13_C15->isAromatic(), true);
    QCOMPARE(C15_C16->isAromatic(), true);
    QCOMPARE(C16_C17->isAromatic(), true);
    QCOMPARE(C16_N20->isAromatic(), true);
    QCOMPARE(C17_C18->isAromatic(), true);
    QCOMPARE(C18_C19->isAromatic(), true);
    QCOMPARE(C19_N20->isAromatic(), true);
    QCOMPARE(C19_C21->isAromatic(), true);
    QCOMPARE(C21_C22->isAromatic(), true);
    QCOMPARE(C22_C23->isAromatic(), true);
    QCOMPARE(C22_N24->isAromatic(), true);
}

/* pyrazole (C3H4N2)
 *
 *      N1
 *    /    \
 *  C5      N2
 *   \\    //
 *   C4 - C3
 */
void RingPerceptionTest::pyrazole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *N2 = molecule.addAtom("N");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *N1_N2 = molecule.addBond(N1, N2, 1);
    chemkit::Bond *N1_C5 = molecule.addBond(N1, C5, 1);
    chemkit::Bond *N2_C3 = molecule.addBond(N2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C3H4N2"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(N2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(N2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(N1_N2->isInRing(), true);
    QCOMPARE(N1_C5->isInRing(), true);
    QCOMPARE(N2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(N1_N2->isAromatic(), true);
    QCOMPARE(N1_C5->isAromatic(), true);
    QCOMPARE(N2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* pyrene (C10H16)
 *
 *           C7  = C8
 *           /       \
 *    C2 = C3        C9 = C14
 *   /       \       /      \
 *  C1       C4  - C10      C15
 *   \\     //      \\      //
 *    C6 - C5       C11 - C16
 *           \       /
 *           C13 = C12
 */
void RingPerceptionTest::pyrene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *C14 = molecule.addAtom("C");
    chemkit::Atom *C15 = molecule.addAtom("C");
    chemkit::Atom *C16 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_C7 = molecule.addBond(C3, C7, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    chemkit::Bond *C4_C10 = molecule.addBond(C4, C10, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    chemkit::Bond *C5_C13 = molecule.addBond(C5, C13, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    chemkit::Bond *C9_C14 = molecule.addBond(C9, C14, 2);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 2);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 1);
    chemkit::Bond *C11_C16 = molecule.addBond(C11, C16, 1);
    chemkit::Bond *C12_C13 = molecule.addBond(C12, C13, 2);
    chemkit::Bond *C14_C15 = molecule.addBond(C14, C15, 1);
    chemkit::Bond *C15_C16 = molecule.addBond(C15, C16, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C16H10"));

    QCOMPARE(molecule.ringCount(), size_t(4));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);
    chemkit::Ring *R3 = C13->smallestRing();
    QCOMPARE(R3->size(), size_t(6));
    QCOMPARE(R3->isAromatic(), true);
    chemkit::Ring *R4 = C14->smallestRing();
    QCOMPARE(R4->size(), size_t(6));
    QCOMPARE(R4->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);
    QCOMPARE(C13->isInRing(), true);
    QCOMPARE(C14->isInRing(), true);
    QCOMPARE(C15->isInRing(), true);
    QCOMPARE(C16->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);
    QCOMPARE(C13->isAromatic(), true);
    QCOMPARE(C14->isAromatic(), true);
    QCOMPARE(C15->isAromatic(), true);
    QCOMPARE(C16->isAromatic(), true);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C7->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C4_C10->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C13->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C9_C14->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);
    QCOMPARE(C11_C16->isInRing(), true);
    QCOMPARE(C12_C13->isInRing(), true);
    QCOMPARE(C14_C15->isInRing(), true);
    QCOMPARE(C15_C16->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C7->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C4_C10->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C5_C13->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C9_C14->isAromatic(), true);
    QCOMPARE(C10_C11->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
    QCOMPARE(C11_C16->isAromatic(), true);
    QCOMPARE(C12_C13->isAromatic(), true);
    QCOMPARE(C14_C15->isAromatic(), true);
    QCOMPARE(C15_C16->isAromatic(), true);
}

/* pyridine (C5H5N)
 *
 *      N1
 *   //    \
 *  C6     C2
 *   |     ||
 *  C5     C3
 *   \\    /
 *      C4
 */
void RingPerceptionTest::pyridine()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Bond *N1_C2 = molecule.addBond(N1, C2, 2);
    chemkit::Bond *N1_C6 = molecule.addBond(N1, C6, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 2);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C5H5N"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);

    QCOMPARE(N1_C2->isInRing(), true);
    QCOMPARE(N1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);

    QCOMPARE(N1_C2->isAromatic(), true);
    QCOMPARE(N1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
}

/* pyrrole (C4H5N)
 *
 *      N1
 *    /    \
 *  C5      C2
 *   \\    //
 *   C4 - C3
 */
void RingPerceptionTest::pyrrole()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *N1_C2 = molecule.addBond(N1, C2, 1);
    chemkit::Bond *N1_C5 = molecule.addBond(N1, C5, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H5N"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(N1_C2->isInRing(), true);
    QCOMPARE(N1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(N1_C2->isAromatic(), true);
    QCOMPARE(N1_C5->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* quinoxaline (C8H6N2)
 *
 *     N1     C7
 *   /   \\ /   \\
 * C6     C2     C8
 * ||      |      |
 * C5     C3     C9
 *   \   // \   //
 *     N4     C10
 */
void RingPerceptionTest::quinoxaline()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *N4 = molecule.addAtom("N");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Bond *N1_C2 = molecule.addBond(N1, C2, 2);
    chemkit::Bond *N1_C6 = molecule.addBond(N1, C6, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_N4 = molecule.addBond(C3, N4, 2);
    chemkit::Bond *C3_C10 = molecule.addBond(C3, C10, 1);
    chemkit::Bond *N4_C5 = molecule.addBond(N4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H6N2"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = N1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), true);

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(N4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(N4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);

    QCOMPARE(N1_C2->isInRing(), true);
    QCOMPARE(N1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_N4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(N4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);

    QCOMPARE(N1_C2->isAromatic(), true);
    QCOMPARE(N1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), true);
    QCOMPARE(C3_N4->isAromatic(), true);
    QCOMPARE(C3_C10->isAromatic(), true);
    QCOMPARE(N4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
}

/* tetralin (C10H12)
 *
 *     C1     C7
 *   /   \\ /    \
 * C6     C2     C8
 * ||      |      |
 * C5     C3     C9
 *   \   // \    /
 *     C4     C10
 */
void RingPerceptionTest::tetralin()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 2);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 1);
    chemkit::Bond *C2_C7 = molecule.addBond(C2, C7, 1);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 2);
    chemkit::Bond *C3_C10 = molecule.addBond(C3, C10, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 2);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 1);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 1);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C10H12"));

    QCOMPARE(molecule.ringCount(), size_t(2));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);
    chemkit::Ring *R2 = C7->smallestRing();
    QCOMPARE(R2->size(), size_t(6));
    QCOMPARE(R2->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);
    QCOMPARE(C9->isAromatic(), false);
    QCOMPARE(C10->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C7->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C10->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C2_C7->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C3_C10->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), true);
    QCOMPARE(C5_C6->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), false);
    QCOMPARE(C8_C9->isAromatic(), false);
    QCOMPARE(C9_C10->isAromatic(), false);
}

/* thiophene (C4H4S)
 *
 *      S1
 *    /    \
 *  C5      C2
 *   \\    //
 *   C4 - C3
 */
void RingPerceptionTest::thiophene()
{
    chemkit::Molecule molecule;
    chemkit::Atom *S1 = molecule.addAtom("S");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Bond *S1_C2 = molecule.addBond(S1, C2, 1);
    chemkit::Bond *S1_C5 = molecule.addBond(S1, C5, 1);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 2);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H4S"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = molecule.rings()[0];
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(S1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);

    QCOMPARE(S1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(C5->isAromatic(), true);

    QCOMPARE(S1_C2->isInRing(), true);
    QCOMPARE(S1_C5->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);

    QCOMPARE(S1_C2->isAromatic(), true);
    QCOMPARE(S1_C5->isAromatic(), true);
    QCOMPARE(C2_C3->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), true);
    QCOMPARE(C4_C5->isAromatic(), true);
}

/* tricyclohexane (C6H8)
 *
 *     /--- C1
 *   C6   /    \
 *   |   C4 -- C2
 *   C5   \    /
 *     \--- C3
 */
void RingPerceptionTest::tricyclohexane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C4 = molecule.addBond(C1, C4);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C2_C4 = molecule.addBond(C2, C4);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C5 = molecule.addBond(C3, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C6H8"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(3));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C3->smallestRing();
    QCOMPARE(R2->size(), size_t(3));
    QCOMPARE(R2->isAromatic(), false);
    chemkit::Ring *R3 = C5->smallestRing();
    QCOMPARE(R3->size(), size_t(5));
    QCOMPARE(R3->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C4->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C4->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C4->isAromatic(), false);
    QCOMPARE(C1_C6->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C2_C4->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C4->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C4->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
}

/* tricyclooctane (C8H12)
 *
 *     C2 - C3
 *    / \   / \
 *   C1  \ /   C4
 *   |    X    |
 *   C8  / \   C5
 *    \ /   \ /
 *     C7 - C6
 */
void RingPerceptionTest::tricyclooctane()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2);
    chemkit::Bond *C1_C8 = molecule.addBond(C1, C8);
    chemkit::Bond *C2_C3 = molecule.addBond(C2, C3);
    chemkit::Bond *C2_C6 = molecule.addBond(C2, C6);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4);
    chemkit::Bond *C3_C7 = molecule.addBond(C3, C7);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6);
    chemkit::Bond *C6_C7 = molecule.addBond(C6, C7);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C8H12"));

    QCOMPARE(molecule.ringCount(), size_t(3));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(5));
    QCOMPARE(R1->isAromatic(), false);
    chemkit::Ring *R2 = C2->smallestRing();
    QCOMPARE(R2->size(), size_t(4));
    QCOMPARE(R2->isAromatic(), false);
    chemkit::Ring *R3 = C4->smallestRing();
    QCOMPARE(R3->size(), size_t(5));
    QCOMPARE(R3->isAromatic(), false);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);

    QCOMPARE(C1->isAromatic(), false);
    QCOMPARE(C2->isAromatic(), false);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), false);
    QCOMPARE(C7->isAromatic(), false);
    QCOMPARE(C8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C8->isInRing(), true);
    QCOMPARE(C2_C3->isInRing(), true);
    QCOMPARE(C2_C6->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C7->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C6_C7->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), false);
    QCOMPARE(C1_C8->isAromatic(), false);
    QCOMPARE(C2_C3->isAromatic(), false);
    QCOMPARE(C2_C6->isAromatic(), false);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C7->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C6_C7->isAromatic(), false);
    QCOMPARE(C7_C8->isAromatic(), false);
}

/* uracil (C4H4N2O2)
 *
 *      O7
 *      ||
 *      C2
 *    /    \
 *  C1      N3
 *  ||      |
 *  C6      C4 = O8
 *    \    /
 *      N5
 */
void RingPerceptionTest::uracil()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *N3 = molecule.addAtom("N");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *N5 = molecule.addAtom("N");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *O7 = molecule.addAtom("O");
    chemkit::Atom *O8 = molecule.addAtom("O");
    chemkit::Bond *C1_C2 = molecule.addBond(C1, C2, 1);
    chemkit::Bond *C1_C6 = molecule.addBond(C1, C6, 2);
    chemkit::Bond *C2_N3 = molecule.addBond(C2, N3, 1);
    chemkit::Bond *C2_O7 = molecule.addBond(C2, O7, 2);
    chemkit::Bond *N3_C4 = molecule.addBond(N3, C4, 1);
    chemkit::Bond *C4_N5 = molecule.addBond(C4, N5, 1);
    chemkit::Bond *C4_O8 = molecule.addBond(C4, O8, 2);
    chemkit::Bond *N5_C6 = molecule.addBond(N5, C6, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C4H4N2O2"));

    QCOMPARE(molecule.ringCount(), size_t(1));
    chemkit::Ring *R1 = C1->smallestRing();
    QCOMPARE(R1->size(), size_t(6));
    QCOMPARE(R1->isAromatic(), true);

    QCOMPARE(C1->isInRing(), true);
    QCOMPARE(C2->isInRing(), true);
    QCOMPARE(N3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(N5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(O7->isInRing(), false);
    QCOMPARE(O8->isInRing(), false);

    QCOMPARE(C1->isAromatic(), true);
    QCOMPARE(C2->isAromatic(), true);
    QCOMPARE(N3->isAromatic(), true);
    QCOMPARE(C4->isAromatic(), true);
    QCOMPARE(N5->isAromatic(), true);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(O7->isAromatic(), false);
    QCOMPARE(O8->isAromatic(), false);

    QCOMPARE(C1_C2->isInRing(), true);
    QCOMPARE(C1_C6->isInRing(), true);
    QCOMPARE(C2_N3->isInRing(), true);
    QCOMPARE(C2_O7->isInRing(), false);
    QCOMPARE(N3_C4->isInRing(), true);
    QCOMPARE(C4_N5->isInRing(), true);
    QCOMPARE(C4_O8->isInRing(), false);
    QCOMPARE(N5_C6->isInRing(), true);

    QCOMPARE(C1_C2->isAromatic(), true);
    QCOMPARE(C1_C6->isAromatic(), true);
    QCOMPARE(C2_N3->isAromatic(), true);
    QCOMPARE(C2_O7->isAromatic(), false);
    QCOMPARE(N3_C4->isAromatic(), true);
    QCOMPARE(C4_N5->isAromatic(), true);
    QCOMPARE(C4_O8->isAromatic(), false);
    QCOMPARE(N5_C6->isAromatic(), true);
}

/* vigtua (C12H8N2)
 *
 *                   N1      C8
 *     -- C5 --\  //    \  /    \\
 *    /    \    C6       C7      C9
 *  C14 -- C4   |        ||      |
 *    \    /    C13      C12     C10
 *     -- C3 --/  \\    /  \    //
 *                   N2      C11
 */
void RingPerceptionTest::vigtua()
{
    chemkit::Molecule molecule;
    chemkit::Atom *N1 = molecule.addAtom("N");
    chemkit::Atom *N2 = molecule.addAtom("N");
    chemkit::Atom *C3 = molecule.addAtom("C");
    chemkit::Atom *C4 = molecule.addAtom("C");
    chemkit::Atom *C5 = molecule.addAtom("C");
    chemkit::Atom *C6 = molecule.addAtom("C");
    chemkit::Atom *C7 = molecule.addAtom("C");
    chemkit::Atom *C8 = molecule.addAtom("C");
    chemkit::Atom *C9 = molecule.addAtom("C");
    chemkit::Atom *C10 = molecule.addAtom("C");
    chemkit::Atom *C11 = molecule.addAtom("C");
    chemkit::Atom *C12 = molecule.addAtom("C");
    chemkit::Atom *C13 = molecule.addAtom("C");
    chemkit::Atom *C14 = molecule.addAtom("C");
    chemkit::Bond *N1_C6 = molecule.addBond(N1, C6, 2);
    chemkit::Bond *N1_C7 = molecule.addBond(N1, C7, 1);
    chemkit::Bond *N2_C12 = molecule.addBond(N2, C12, 1);
    chemkit::Bond *N2_C13 = molecule.addBond(N2, C13, 2);
    chemkit::Bond *C3_C4 = molecule.addBond(C3, C4, 1);
    chemkit::Bond *C3_C13 = molecule.addBond(C3, C13, 1);
    chemkit::Bond *C3_C14 = molecule.addBond(C3, C14, 1);
    chemkit::Bond *C4_C5 = molecule.addBond(C4, C5, 1);
    chemkit::Bond *C4_C14 = molecule.addBond(C4, C14, 1);
    chemkit::Bond *C5_C6 = molecule.addBond(C5, C6, 1);
    chemkit::Bond *C5_C14 = molecule.addBond(C5, C14, 1);
    chemkit::Bond *C6_C13 = molecule.addBond(C6, C13, 1);
    chemkit::Bond *C7_C8 = molecule.addBond(C7, C8, 1);
    chemkit::Bond *C7_C12 = molecule.addBond(C7, C12, 2);
    chemkit::Bond *C8_C9 = molecule.addBond(C8, C9, 2);
    chemkit::Bond *C9_C10 = molecule.addBond(C9, C10, 1);
    chemkit::Bond *C10_C11 = molecule.addBond(C10, C11, 2);
    chemkit::Bond *C11_C12 = molecule.addBond(C11, C12, 1);
    addHydrogens(&molecule);
    QCOMPARE(molecule.formula(), std::string("C12H8N2"));

    QCOMPARE(molecule.ringCount(), size_t(5));
    chemkit::Ring *R1 = C5->smallestRing();
    QCOMPARE(R1->size(), size_t(3));
    chemkit::Ring *R2 = C3->smallestRing();
    QCOMPARE(R2->size(), size_t(3));
    chemkit::Ring *R3 = C6->smallestRing();
    QCOMPARE(R3->size(), size_t(5));
    chemkit::Ring *R4 = N1->smallestRing();
    QCOMPARE(R4->size(), size_t(6));
    chemkit::Ring *R5 = C8->smallestRing();
    QCOMPARE(R5->size(), size_t(6));

    QCOMPARE(N1->isInRing(), true);
    QCOMPARE(N2->isInRing(), true);
    QCOMPARE(C3->isInRing(), true);
    QCOMPARE(C4->isInRing(), true);
    QCOMPARE(C5->isInRing(), true);
    QCOMPARE(C6->isInRing(), true);
    QCOMPARE(C7->isInRing(), true);
    QCOMPARE(C8->isInRing(), true);
    QCOMPARE(C9->isInRing(), true);
    QCOMPARE(C10->isInRing(), true);
    QCOMPARE(C11->isInRing(), true);
    QCOMPARE(C12->isInRing(), true);
    QCOMPARE(C13->isInRing(), true);
    QCOMPARE(C14->isInRing(), true);

    QCOMPARE(N1->isAromatic(), true);
    QCOMPARE(N2->isAromatic(), true);
    QCOMPARE(C3->isAromatic(), false);
    QCOMPARE(C4->isAromatic(), false);
    QCOMPARE(C5->isAromatic(), false);
    QCOMPARE(C6->isAromatic(), true);
    QCOMPARE(C7->isAromatic(), true);
    QCOMPARE(C8->isAromatic(), true);
    QCOMPARE(C9->isAromatic(), true);
    QCOMPARE(C10->isAromatic(), true);
    QCOMPARE(C11->isAromatic(), true);
    QCOMPARE(C12->isAromatic(), true);
    QCOMPARE(C13->isAromatic(), true);
    QCOMPARE(C14->isAromatic(), false);

    QCOMPARE(N1_C6->isInRing(), true);
    QCOMPARE(N1_C7->isInRing(), true);
    QCOMPARE(N2_C12->isInRing(), true);
    QCOMPARE(N2_C13->isInRing(), true);
    QCOMPARE(C3_C4->isInRing(), true);
    QCOMPARE(C3_C13->isInRing(), true);
    QCOMPARE(C3_C14->isInRing(), true);
    QCOMPARE(C4_C5->isInRing(), true);
    QCOMPARE(C4_C14->isInRing(), true);
    QCOMPARE(C5_C6->isInRing(), true);
    QCOMPARE(C5_C14->isInRing(), true);
    QCOMPARE(C6_C13->isInRing(), true);
    QCOMPARE(C7_C8->isInRing(), true);
    QCOMPARE(C7_C12->isInRing(), true);
    QCOMPARE(C8_C9->isInRing(), true);
    QCOMPARE(C9_C10->isInRing(), true);
    QCOMPARE(C10_C11->isInRing(), true);
    QCOMPARE(C11_C12->isInRing(), true);

    QCOMPARE(N1_C6->isAromatic(), true);
    QCOMPARE(N1_C7->isAromatic(), true);
    QCOMPARE(N2_C12->isAromatic(), true);
    QCOMPARE(N2_C13->isAromatic(), true);
    QCOMPARE(C3_C4->isAromatic(), false);
    QCOMPARE(C3_C13->isAromatic(), false);
    QCOMPARE(C3_C14->isAromatic(), false);
    QCOMPARE(C4_C5->isAromatic(), false);
    QCOMPARE(C4_C14->isAromatic(), false);
    QCOMPARE(C5_C6->isAromatic(), false);
    QCOMPARE(C5_C14->isAromatic(), false);
    QCOMPARE(C6_C13->isAromatic(), true);
    QCOMPARE(C7_C8->isAromatic(), true);
    QCOMPARE(C7_C12->isAromatic(), true);
    QCOMPARE(C8_C9->isAromatic(), true);
    QCOMPARE(C9_C10->isAromatic(), true);
    QCOMPARE(C10_C11->isAromatic(), true);
    QCOMPARE(C11_C12->isAromatic(), true);
}

QTEST_APPLESS_MAIN(RingPerceptionTest)
