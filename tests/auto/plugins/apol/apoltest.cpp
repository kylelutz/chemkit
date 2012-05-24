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

#include "apoltest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void ApolTest::initTestCase()
{
    // verify that the apol plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "apol") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "bpol") == 1);
}

void ApolTest::apol_data()
{
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("apol");
    QTest::addColumn<double>("bpol");

    QTest::newRow("methane") << "C" << "CH4" << 4.427172 << 4.372828;
    QTest::newRow("ammonia") << "N" << "H3N" << 3.10038 << 1.29962;
    QTest::newRow("ethanol") << "CCO" << "C2H6O" << 8.322758 << 6.559242;
}

void ApolTest::apol()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);
    QFETCH(double, apol);
    QFETCH(double, bpol);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());
    QCOMPARE(qRound(molecule.descriptor("apol").toDouble() * 100), qRound(apol * 100));
    QCOMPARE(qRound(molecule.descriptor("bpol").toDouble() * 100), qRound(bpol * 100));
}

QTEST_APPLESS_MAIN(ApolTest)
