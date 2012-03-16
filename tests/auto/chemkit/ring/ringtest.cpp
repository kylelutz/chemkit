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
    QCOMPARE(benzene->ringCount(), size_t(1));
    benzeneRing = benzene->rings()[0];

    // furan
    furan = inchi->read("InChI=1/C4H4O/c1-2-4-5-3-1/h1-4H");
    QVERIFY(furan != 0);
    QCOMPARE(furan->ringCount(), size_t(1));
    furanRing = furan->rings()[0];

    // cyclohexane
    cyclohexane = inchi->read("InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2");
    QVERIFY(cyclohexane != 0);
    QCOMPARE(cyclohexane->ringCount(), size_t(1));
    cyclohexaneRing = cyclohexane->rings()[0];

    // cyclopropane
    cyclopropane = inchi->read("InChI=1/C3H6/c1-2-3-1/h1-3H2");
    QVERIFY(cyclopropane != 0);
    QCOMPARE(cyclopropane->ringCount(), size_t(1));
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
    QCOMPARE(size_t(benzeneRing->atoms().size()), size_t(6));
    foreach(chemkit::Atom *atom, benzeneRing->atoms()){
        QVERIFY(atom->molecule() == benzene);
        QVERIFY(benzene->contains(atom));

        QVERIFY(atom->is(chemkit::Atom::Carbon));
    }
}

void RingTest::atomCount()
{
    QCOMPARE(benzeneRing->atomCount(), size_t(6));
    QCOMPARE(furanRing->atomCount(), size_t(5));
    QCOMPARE(cyclohexaneRing->atomCount(), size_t(6));
    QCOMPARE(cyclopropaneRing->atomCount(), size_t(3));

    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Carbon), size_t(6));
    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Hydrogen), size_t(0));
    QCOMPARE(benzeneRing->atomCount(chemkit::Atom::Oxygen), size_t(0));
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Carbon), size_t(4));
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Hydrogen), size_t(0));
    QCOMPARE(furanRing->atomCount(chemkit::Atom::Oxygen), size_t(1));
}

void RingTest::size()
{
    QCOMPARE(benzeneRing->size(), size_t(6));
    QCOMPARE(furanRing->size(), size_t(5));
    QCOMPARE(cyclohexaneRing->size(), size_t(6));
    QCOMPARE(cyclopropaneRing->size(), size_t(3));
}

void RingTest::bonds()
{
    foreach(chemkit::Bond *bond, benzeneRing->bonds()){
        QVERIFY(bond->molecule() == benzene);
        QVERIFY(benzene->contains(bond));

        QVERIFY(bond->atom1()->is(chemkit::Atom::Carbon));
        QVERIFY(bond->atom2()->is(chemkit::Atom::Carbon));
    }
}

void RingTest::bondCount()
{
    QCOMPARE(benzeneRing->bondCount(), size_t(6));
    QCOMPARE(furanRing->bondCount(), size_t(5));
    QCOMPARE(cyclohexaneRing->bondCount(), size_t(6));
    QCOMPARE(cyclopropaneRing->bondCount(), size_t(3));
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
            QCOMPARE(furanRing->position(atom), size_t(1));
            QCOMPARE(furanRing->position(atom, furanOxygen), size_t(1));
        }
        else if(atom == furanOxygen){
            QCOMPARE(furanRing->position(atom), size_t(0));
            QCOMPARE(furanRing->position(atom, furanOxygen), size_t(0));
        }
        else{
            QCOMPARE(furanRing->position(atom), size_t(2));
            QCOMPARE(furanRing->position(atom, furanOxygen), size_t(2));
        }
    }

    // atom from another molecule
    QCOMPARE(furanRing->position(benzene->atoms()[0]), size_t(0));
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
    QCOMPARE(benzeneRing->heteroatomCount(), size_t(0));
    QCOMPARE(furanRing->heteroatomCount(), size_t(1));
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
