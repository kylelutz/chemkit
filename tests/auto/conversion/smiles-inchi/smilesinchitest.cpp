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

#include "smilesinchitest.h"

#include <algorithm>

#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

void SmilesInchiTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::LineFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "smiles") != formats.end());
    QVERIFY(std::find(formats.begin(), formats.end(), "inchi") != formats.end());
}

void SmilesInchiTest::alanine()
{
    // D-alanine
    chemkit::Molecule molecule("O=C(O)[C@H](N)C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1"));

    // L-alanine
    molecule = chemkit::Molecule("C[C@@H](C(=O)O)N", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"));
}

void SmilesInchiTest::butene()
{
    // cis-butene
    chemkit::Molecule molecule("C/C=C\\C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3-"));

    // trans-butene
    molecule = chemkit::Molecule("C/C=C/C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+"));
}

void SmilesInchiTest::carvone()
{
    // S-Carvone
    chemkit::Molecule molecule("O=C1\\C(=C/C[C@H](/C(=C)C)C1)C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C10H14O/c1-7(2)9-5-4-8(3)10(11)6-9/h4,9H,1,5-6H2,2-3H3/t9-/m0/s1"));

    // R-Carvone
    molecule = chemkit::Molecule("O=C1\\C(=C/C[C@@H](/C(=C)C)C1)C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C10H14O/c1-7(2)9-5-4-8(3)10(11)6-9/h4,9H,1,5-6H2,2-3H3/t9-/m1/s1"));
}

void SmilesInchiTest::ethanol()
{
    chemkit::Molecule molecule("CCO", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"));
}

void SmilesInchiTest::serine()
{
    // D-serine
    chemkit::Molecule molecule("O=C(O)[C@H](N)CO", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)/t2-/m1/s1"));

    // L-serine
    molecule = chemkit::Molecule("C([C@@H](C(=O)O)N)O", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)/t2-/m0/s1"));
}

void SmilesInchiTest::threonine()
{
    // L-threonine
    chemkit::Molecule molecule("C[C@H]([C@@H](C(=O)O)N)O", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/t2-,3+/m1/s1"));

    // D-threonine
    molecule = chemkit::Molecule("O=C(O)[C@H](N)[C@@H](O)C", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C4H9NO3/c1-2(6)3(5)4(7)8/h2-3,6H,5H2,1H3,(H,7,8)/t2-,3+/m0/s1"));
}

void SmilesInchiTest::tyrosine()
{
    // L-tyrosine
    chemkit::Molecule molecule("c1cc(ccc1C[C@@H](C(=O)O)N)O", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4-2-6/h1-4,8,11H,5,10H2,(H,12,13)/t8-/m0/s1"));

    // D-tyrosine
    molecule = chemkit::Molecule("O=C(O)[C@H](N)Cc1ccc(O)cc1", "smiles");
    QCOMPARE(molecule.formula("inchi"),
             std::string("InChI=1S/C9H11NO3/c10-8(9(12)13)5-6-1-3-7(11)4-2-6/h1-4,8,11H,5,10H2,(H,12,13)/t8-/m1/s1"));
}

QTEST_APPLESS_MAIN(SmilesInchiTest)
