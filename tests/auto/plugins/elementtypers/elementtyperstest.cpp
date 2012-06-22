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

#include "elementtyperstest.h"

#include <boost/scoped_ptr.hpp>
#include <boost/range/algorithm.hpp>

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>

void ElementTypersTest::initTestCase()
{
    // verify that the elementtypers plugin registered itself correctly
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "element-name") == 1);
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "atomic-number") == 1);
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "element-symbol") == 1);
}

void ElementTypersTest::test_data()
{
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");

    QTest::newRow("ethanol") << "CCO" << "C2H6O";
    QTest::newRow("guanine") << "NC1=Nc2[nH]cnc2C(=O)N1" << "C5H5N5O";
}

void ElementTypersTest::test()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());

    boost::scoped_ptr<chemkit::AtomTyper> elementNameTyper(chemkit::AtomTyper::create("element-name"));
    QVERIFY(elementNameTyper != 0);
    QCOMPARE(elementNameTyper->name(), std::string("element-name"));
    elementNameTyper->setMolecule(&molecule);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        QCOMPARE(elementNameTyper->type(atom), atom->name());
    }

    boost::scoped_ptr<chemkit::AtomTyper> atomicNumberTyper(chemkit::AtomTyper::create("atomic-number"));
    QVERIFY(atomicNumberTyper != 0);
    QCOMPARE(atomicNumberTyper->name(), std::string("atomic-number"));
    atomicNumberTyper->setMolecule(&molecule);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        QCOMPARE(boost::lexical_cast<int>(atomicNumberTyper->type(atom)), int(atom->atomicNumber()));
    }

    boost::scoped_ptr<chemkit::AtomTyper> elementSymbolTyper(chemkit::AtomTyper::create("element-symbol"));
    QVERIFY(elementSymbolTyper != 0);
    QCOMPARE(elementSymbolTyper->name(), std::string("element-symbol"));
    elementSymbolTyper->setMolecule(&molecule);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        QCOMPARE(elementSymbolTyper->type(atom), atom->symbol());
    }
}

QTEST_APPLESS_MAIN(ElementTypersTest)
