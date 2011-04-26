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

#include "bondpredictortest.h"

#include <chemkit/molecule.h>
#include <chemkit/bondpredictor.h>

void BondPredictorTest::predictBonds()
{
    // create di-hydrogen molecule
    chemkit::Molecule hydrogenMolecule;
    chemkit::Atom *h1 = hydrogenMolecule.addAtom(1);
    chemkit::Atom *h2 = hydrogenMolecule.addAtom(1);

    // set atoms 1 angstrom apart
    h1->setPosition(0, 0, 0);
    h2->setPosition(1, 0, 0);

    // predict bonds and verify the hydrogens are bonded
    chemkit::BondPredictor::predictBonds(&hydrogenMolecule);
    QCOMPARE(hydrogenMolecule.bondCount(), 1);
    QCOMPARE(h1->isBondedTo(h2), true);

    // remove bond, move 100 angstroms apart
    hydrogenMolecule.removeBond(h1, h2);
    h2->setPosition(100, 0, 0);

    // predict bonds and verify the hydrogens are not bonded
    chemkit::BondPredictor::predictBonds(&hydrogenMolecule);
    QCOMPARE(hydrogenMolecule.bondCount(), 0);
    QCOMPARE(h1->isBondedTo(h2), false);
}

QTEST_APPLESS_MAIN(BondPredictorTest)
