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

#include "ringtest.h"

#include <chemkit/chemkit.h>
#include <chemkit/ring.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void RingTest::initTestCase()
{
    chemkit::LineFormat *inchi = chemkit::LineFormat::create("inchi");

    // benzene
    benzene = inchi->read("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H");
    QVERIFY(benzene != 0);
    QCOMPARE(benzene->ringCount(), 1);
    benzeneRing = benzene->rings()[0];

    // furan
    furan = inchi->read("InChI=1/C4H4O/c1-2-4-5-3-1/h1-4H");
    QVERIFY(furan != 0);
    QCOMPARE(furan->ringCount(), 1);
    furanRing = furan->rings()[0];

    // cyclohexane
    cyclohexane = inchi->read("InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2");
    QVERIFY(cyclohexane != 0);
    QCOMPARE(cyclohexane->ringCount(), 1);
    cyclohexaneRing = cyclohexane->rings()[0];

    // cyclopropane
    cyclopropane = inchi->read("InChI=1/C3H6/c1-2-3-1/h1-3H2");
    QVERIFY(cyclopropane != 0);
    QCOMPARE(cyclopropane->ringCount(), 1);
    cyclopropaneRing = cyclopropane->rings()[0];

    delete inchi;
}

void RingTest::molecule()
{
    QCOMPARE(benzeneRing->molecule(), benzene);
    QCOMPARE(furanRing->molecule(), furan);
}

void RingTest::atoms()
{
    QCOMPARE(benzeneRing->atoms().size(), 6UL);
    foreach(chemkit::Atom *atom, benzeneRing->atoms()){
        QVERIFY(atom->molecule() == benzene);
        QVERIFY(benzene->contains(atom));

        QVERIFY(atom->is(chemkit::Atom::Carbon));
    }
}

void RingTest::atomCount()
{
    QCOMPARE(benzeneRing->atomCount(), 6);
    QCOMPARE(furanRing->atomCount(), 5);
    QCOMPARE(cyclohexaneRing->atomCount(), 6);
    QCOMPARE(cyclopropaneRing->atomCount(), 3);

    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Carbon), 6);
    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Hydrogen), 0);
    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Oxygen), 0);
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Carbon), 4);
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Hydrogen), 0);
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Oxygen), 1);
}

void RingTest::size()
{
    QCOMPARE(benzeneRing->size(), 6);
    QCOMPARE(furanRing->size(), 5);
    QCOMPARE(cyclohexaneRing->size(), 6);
    QCOMPARE(cyclopropaneRing->size(), 3);
}

void RingTest::bonds()
{
    QCOMPARE(benzeneRing->bonds().size(), 6);
    foreach(chemkit::Bond *bond, benzeneRing->bonds()){
        QVERIFY(bond->molecule() == benzene);
        QVERIFY(benzene->contains(bond));

        QVERIFY(bond->atom1()->is(chemkit::Atom::Carbon));
        QVERIFY(bond->atom2()->is(chemkit::Atom::Carbon));
    }
}

void RingTest::bondCount()
{
    QCOMPARE(benzeneRing->bondCount(), 6);
    QCOMPARE(furanRing->bondCount(), 5);
    QCOMPARE(cyclohexaneRing->bondCount(), 6);
    QCOMPARE(cyclopropaneRing->bondCount(), 3);
}

void RingTest::root()
{
    QVERIFY(benzeneRing->root()->is(chemkit::Atom::Carbon));
    QVERIFY(furanRing->root()->is(chemkit::Atom::Oxygen));
}

void RingTest::position()
{
    chemkit::Atom *furanOxygen = furanRing->root();
    QVERIFY(furanOxygen->is(chemkit::Atom::Oxygen));
    foreach(chemkit::Atom *atom, furanRing->atoms()){
        if(atom->isBondedTo(furanOxygen)){
            QCOMPARE(furanRing->position(atom), 1);
            QCOMPARE(furanRing->position(atom, furanOxygen), 1);
        }
        else if(atom == furanOxygen){
            QCOMPARE(furanRing->position(atom), 0);
            QCOMPARE(furanRing->position(atom, furanOxygen), 0);
        }
        else{
            QCOMPARE(furanRing->position(atom), 2);
            QCOMPARE(furanRing->position(atom, furanOxygen), 2);
        }
    }

    // atom from another molecule
    QCOMPARE(furanRing->position(benzene->atoms()[0]), 0);
}

void RingTest::contains()
{
    foreach(chemkit::Atom *atom, benzene->atoms()){
        if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(benzeneRing->contains(atom), true);
        }
        else{
            QCOMPARE(benzeneRing->contains(atom), false);
        }
    }
    foreach(chemkit::Bond *bond, benzene->bonds()){
        if(bond->atom1()->is(chemkit::Atom::Carbon) && bond->atom2()->is(chemkit::Atom::Carbon)){
            QCOMPARE(benzeneRing->contains(bond), true);
        }
        else{
            QCOMPARE(benzeneRing->contains(bond), false);
        }
    }
    QCOMPARE(benzeneRing->contains(chemkit::Atom::Carbon), true);
    QCOMPARE(benzeneRing->contains(chemkit::Atom::Hydrogen), false);
    QCOMPARE(benzeneRing->contains(chemkit::Atom::Oxygen), false);

    QCOMPARE(furanRing->contains(chemkit::Atom::Carbon), true);
    QCOMPARE(furanRing->contains(chemkit::Atom::Hydrogen), false);
    QCOMPARE(furanRing->contains(chemkit::Atom::Oxygen), true);
}

void RingTest::heteroatomCount()
{
    QCOMPARE(benzeneRing->heteroatomCount(), 0);
    QCOMPARE(furanRing->heteroatomCount(), 1);
}

void RingTest::isHeterocycle()
{
    QCOMPARE(benzeneRing->isHeterocycle(), false);
    QCOMPARE(furanRing->isHeterocycle(), true);
}

void RingTest::isAromatic()
{
    QCOMPARE(benzeneRing->isAromatic(), true);
    QCOMPARE(furanRing->isAromatic(), true);
    QCOMPARE(cyclohexaneRing->isAromatic(), false);
    QCOMPARE(cyclopropaneRing->isAromatic(), false);
}

void RingTest::cleanupTestCase()
{
    delete benzeneRing->molecule();
    delete furanRing->molecule();
    delete cyclohexaneRing->molecule();
    delete cyclopropaneRing->molecule();
}

QTEST_APPLESS_MAIN(RingTest)
