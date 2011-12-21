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

#include "moietytest.h"

#include <chemkit/moiety.h>
#include <chemkit/molecule.h>
#include <vector>

void MoietyTest::size()
{
    chemkit::Molecule molecule;
    chemkit::Atom* atomH = molecule.addAtom("H");
    chemkit::Atom* atomC = molecule.addAtom("C");
    std::vector<chemkit::Atom *> moietyAtoms;
    moietyAtoms.push_back(atomH);
    moietyAtoms.push_back(atomC);
    chemkit::Moiety moiety(moietyAtoms);
    QCOMPARE(moiety.size(), size_t(2));

    QCOMPARE(chemkit::Moiety().size(), size_t(0));
}

void MoietyTest::isEmpty()
{
    chemkit::Molecule* molecule = new chemkit::Molecule;
    chemkit::Atom* atomH = molecule->addAtom("H");
    chemkit::Atom* atomC = molecule->addAtom("C");
    std::vector<chemkit::Atom*> moietyAtoms;
    moietyAtoms.push_back(atomH);
    moietyAtoms.push_back(atomC);
    chemkit::Moiety moiety(moietyAtoms);
    QCOMPARE(moiety.isEmpty(), false);

    QCOMPARE(chemkit::Moiety().isEmpty(),true);

    delete molecule;
}

void MoietyTest::atomCount()
{
    chemkit::Molecule* molecule = new chemkit::Molecule;
    chemkit::Atom* atomH = molecule->addAtom("H");
    chemkit::Atom* atomC = molecule->addAtom("C");
    std::vector<chemkit::Atom*> moietyAtoms;
    moietyAtoms.push_back(atomH);
    moietyAtoms.push_back(atomC);
    chemkit::Moiety moiety(moietyAtoms);
    QCOMPARE(moiety.atomCount(), size_t(2));

    QCOMPARE(chemkit::Moiety().atomCount(), size_t(0));

    delete molecule;
}

void MoietyTest::molecule()
{
    chemkit::Molecule* molecule = new chemkit::Molecule;
    chemkit::Atom* atomH = molecule->addAtom("H");
    std::vector<chemkit::Atom*> moietyAtoms;
    moietyAtoms.push_back(atomH);
    chemkit::Moiety moiety(moietyAtoms);
    QVERIFY(moiety.molecule() == molecule);

    chemkit::Moiety emptyMoiety;
    QVERIFY(emptyMoiety.molecule() == 0);

    delete molecule;
}

QTEST_APPLESS_MAIN(MoietyTest)
