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

#include "similaritydescriptortest.h"

#include <chemkit/molecule.h>
#include <chemkit/similaritydescriptor.h>

void SimilarityDescriptorTest::name()
{
    chemkit::SimilarityDescriptor descriptor;
    QCOMPARE(descriptor.name(), std::string("similarity"));
}

void SimilarityDescriptorTest::molecule()
{
    chemkit::SimilarityDescriptor descriptor;
    QVERIFY(descriptor.molecule() == 0);

    chemkit::Molecule molecule;
    descriptor.setMolecule(&molecule);
    QVERIFY(descriptor.molecule() == &molecule);

    descriptor.setMolecule(0);
    QVERIFY(descriptor.molecule() == 0);
}

void SimilarityDescriptorTest::fingerprint()
{
    chemkit::SimilarityDescriptor descriptor;
    QCOMPARE(descriptor.fingerprint(), std::string("fp2"));

    descriptor.setFingerprint(std::string());
    QVERIFY(descriptor.fingerprint() == std::string());

    descriptor.setFingerprint("fp2");
    QVERIFY(descriptor.fingerprint() == std::string("fp2"));
}

void SimilarityDescriptorTest::value()
{
    chemkit::Molecule ethanol("CCO", "smiles");
    QCOMPARE(ethanol.formula(), std::string("C2H6O"));

    chemkit::SimilarityDescriptor descriptor;
    QCOMPARE(qRound(descriptor.value(&ethanol).toDouble()), 0);

    descriptor.setMolecule(&ethanol);
    QCOMPARE(qRound(descriptor.value(&ethanol).toDouble()), 1);

    chemkit::Molecule methanol("CO", "smiles");
    QCOMPARE(qRound(descriptor.value(&methanol).toDouble() * 100), 33);
}

QTEST_APPLESS_MAIN(SimilarityDescriptorTest)
