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
