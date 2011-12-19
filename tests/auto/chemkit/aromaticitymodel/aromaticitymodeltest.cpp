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

#include "aromaticitymodeltest.h"

#include <chemkit/ring.h>
#include <chemkit/molecule.h>

#include "mockaromaticitymodel.h"

void AromaticityModelTest::name()
{
    MockAromaticityModel model;
    QCOMPARE(model.name(), std::string("mock"));
}

void AromaticityModelTest::setMolecule()
{
    MockAromaticityModel model;
    QVERIFY(model.molecule() == 0);

    chemkit::Molecule molecule;
    model.setMolecule(&molecule);
    QVERIFY(model.molecule() == &molecule);

    model.setMolecule(0);
    QVERIFY(model.molecule() == 0);
}

void AromaticityModelTest::isAromatic()
{
    MockAromaticityModel model;

    chemkit::Molecule furan("c1ccoc1", "smiles");
    QCOMPARE(furan.formula(), std::string("C4H4O"));
    QCOMPARE(furan.ringCount(), 1);
    chemkit::Ring *ring = furan.ring(0);
    QCOMPARE(ring->size(), 5);
    QCOMPARE(model.isAromatic(ring), true);
    QCOMPARE(model.isAromatic(ring->atom(0)), true);
    QCOMPARE(model.isAromatic(ring->bond(0)), true);

    chemkit::Molecule cyclohexane("C1CCCCC1", "smiles");
    QCOMPARE(cyclohexane.formula(), std::string("C6H12"));
    QCOMPARE(cyclohexane.ringCount(), 1);
    ring = cyclohexane.ring(0);
    QCOMPARE(ring->size(), 6);
    QCOMPARE(model.isAromatic(ring), false);
    QCOMPARE(model.isAromatic(ring->atom(0)), false);
    QCOMPARE(model.isAromatic(ring->bond(0)), false);
}

QTEST_APPLESS_MAIN(AromaticityModelTest)
