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

#include "randicindextest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void RandicIndexTest::initTestCase()
{
    // verify that the randicindex plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "randic-index") == 1);
}

void RandicIndexTest::ethane()
{
    chemkit::Molecule ethane("CC", "smiles");
    QCOMPARE(ethane.formula(), std::string("C2H6"));

    // index = 1.0
    QCOMPARE(ethane.descriptor("randic-index").toInt(), 1);
}

void RandicIndexTest::isobutane()
{
    chemkit::Molecule isobutane("CC(C)C", "smiles");
    QCOMPARE(isobutane.formula(), std::string("C4H10"));

    // index = 1.7321
    QCOMPARE(qRound(isobutane.descriptor("randic-index").toDouble()), 2);
}

void RandicIndexTest::dimethylpropane()
{
    chemkit::Molecule dimethylpropane("CC(C)(C)C", "smiles");
    QCOMPARE(dimethylpropane.formula(), std::string("C5H12"));

    // index = 2.0
    QCOMPARE(qRound(dimethylpropane.descriptor("randic-index").toDouble()), 2);
}

void RandicIndexTest::octane()
{
    chemkit::Molecule octane("CCCCCCCC", "smiles");
    QCOMPARE(octane.formula(), std::string("C8H18"));

    // index = 3.9142
    QCOMPARE(qRound(octane.descriptor("randic-index").toDouble()), 4);
}

QTEST_APPLESS_MAIN(RandicIndexTest)
