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

#include "hydrogenbondstest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void HydrogenBondsTest::initTestCase()
{
    // verify that the hydrogenbonds plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "hydrogen-bond-donors") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "hydrogen-bond-acceptors") == 1);
}

void HydrogenBondsTest::ethanol()
{
    chemkit::Molecule ethanol("CCO", "smiles");
    QCOMPARE(ethanol.formula(), std::string("C2H6O"));

    QCOMPARE(ethanol.descriptor("hydrogen-bond-donors").toInt(), 1);
    QCOMPARE(ethanol.descriptor("hydrogen-bond-acceptors").toInt(), 1);
}

void HydrogenBondsTest::guanine()
{
    chemkit::Molecule guanine("c1[nH]c2c(n1)c(=O)[nH]c(n2)N", "smiles");
    QCOMPARE(guanine.formula(), std::string("C5H5N5O"));

    QCOMPARE(guanine.descriptor("hydrogen-bond-donors").toInt(), 4);
    QCOMPARE(guanine.descriptor("hydrogen-bond-acceptors").toInt(), 6);
}

QTEST_APPLESS_MAIN(HydrogenBondsTest)
