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

#include "conformertest.h"

#include <chemkit/molecule.h>
#include <chemkit/conformer.h>

void ConformerTest::basic()
{
    chemkit::Molecule molecule;
    chemkit::Conformer *conformer = molecule.conformer();
    QVERIFY(conformer->molecule() == &molecule);

    conformer = molecule.addConformer();
    QVERIFY(conformer->molecule() == &molecule);
}

void ConformerTest::atomPosition()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("C");
    atom->setPosition(1, 2, 3);
    QCOMPARE(atom->position(), chemkit::Point3(1, 2, 3));

    chemkit::Conformer *conformer = molecule.addConformer();
    QVERIFY(conformer != 0);

    conformer->setPosition(atom, chemkit::Point3(3, 2, 1));
    QCOMPARE(conformer->position(atom), chemkit::Point3(3, 2, 1));
    QCOMPARE(atom->position(), chemkit::Point3(1, 2, 3));

    chemkit::Conformer *originalConformer = molecule.conformer();
    QVERIFY(originalConformer != 0);

    molecule.setConformer(conformer);
    QCOMPARE(atom->position(), chemkit::Point3(3, 2, 1));
    QCOMPARE(originalConformer->position(atom), chemkit::Point3(1, 2, 3));
}

QTEST_APPLESS_MAIN(ConformerTest)
