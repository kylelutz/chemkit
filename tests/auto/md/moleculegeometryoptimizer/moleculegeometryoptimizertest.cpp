/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "moleculegeometryoptimizertest.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculegeometryoptimizer.h>

void MoleculeGeometryOptimizerTest::molecule()
{
    chemkit::MoleculeGeometryOptimizer optimizer;
    QVERIFY(optimizer.molecule() == 0);

    chemkit::Molecule molecule;
    optimizer.setMolecule(&molecule);
    QVERIFY(optimizer.molecule() == &molecule);

    optimizer.setMolecule(0);
    QVERIFY(optimizer.molecule() == 0);
}

void MoleculeGeometryOptimizerTest::water()
{
    // build water molecule
    chemkit::Molecule molecule;
    chemkit::Atom *O1 = molecule.addAtom("O");
    chemkit::Atom *H2 = molecule.addAtom("H");
    chemkit::Atom *H3 = molecule.addAtom("H");
    molecule.addBond(O1, H2);
    molecule.addBond(O1, H3);
    QCOMPARE(molecule.formula(), std::string("H2O"));

    O1->setPosition(0, 0, 0);
    H2->setPosition(0, 1, 0);
    H3->setPosition(1, 0, 0);
    QCOMPARE(qRound(molecule.bondAngle(H2, O1, H3)), 90);

    // setup geometry optimizer
    chemkit::MoleculeGeometryOptimizer optimizer(&molecule);
    QVERIFY(optimizer.molecule() == &molecule);

    // optimize
    bool ok = optimizer.optimize();
    if(!ok)
        qDebug() << optimizer.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(qRound(molecule.bondAngle(H2, O1, H3)), 104);
}

QTEST_APPLESS_MAIN(MoleculeGeometryOptimizerTest)
