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

#include "structuresimilaritydescriptortest.h"

#include <boost/make_shared.hpp>

#include <chemkit/molecule.h>
#include <chemkit/structuresimilaritydescriptor.h>

void StructureSimilarityDescriptorTest::name()
{
    chemkit::StructureSimilarityDescriptor descriptor;
    QCOMPARE(descriptor.name(), std::string("structure-similarity"));
}

void StructureSimilarityDescriptorTest::molecule()
{
    chemkit::StructureSimilarityDescriptor descriptor;
    QVERIFY(descriptor.molecule() == boost::shared_ptr<chemkit::Molecule>());

    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    descriptor.setMolecule(molecule);
    QVERIFY(descriptor.molecule() == molecule);

    descriptor.setMolecule(boost::shared_ptr<chemkit::Molecule>());
    QVERIFY(descriptor.molecule() == boost::shared_ptr<chemkit::Molecule>());
}

void StructureSimilarityDescriptorTest::value()
{
    boost::shared_ptr<chemkit::Molecule> methanol =
        boost::make_shared<chemkit::Molecule>("CO", "smiles");
    QCOMPARE(methanol->formula(), std::string("CH4O"));

    boost::shared_ptr<chemkit::Molecule> ethanol =
        boost::make_shared<chemkit::Molecule>("CCO", "smiles");
    QCOMPARE(ethanol->formula(), std::string("C2H6O"));

    boost::shared_ptr<chemkit::Molecule> propanol =
        boost::make_shared<chemkit::Molecule>("CCCO", "smiles");
    QCOMPARE(propanol->formula(), std::string("C3H8O"));

    chemkit::StructureSimilarityDescriptor descriptor;

    // methanol -> methanol
    descriptor.setMolecule(methanol);
    QCOMPARE(descriptor.value(methanol.get()).toInt(), 1);

    // methanol -> ethanol
    QCOMPARE(qRound(descriptor.value(ethanol.get()).toDouble() * 100), 67);

    // methanol -> propanol
    QCOMPARE(qRound(descriptor.value(propanol.get()).toDouble() * 100), 50);

    // ethanol -> ethanol
    descriptor.setMolecule(ethanol);
    QCOMPARE(descriptor.value(ethanol.get()).toInt(), 1);

    // ethanol -> methanol
    QCOMPARE(qRound(descriptor.value(methanol.get()).toDouble() * 100), 67);

    // ethanol -> propanol
    QCOMPARE(qRound(descriptor.value(propanol.get()).toDouble() * 100), 75);

    // propanol -> propanol
    descriptor.setMolecule(propanol);
    QCOMPARE(descriptor.value(propanol.get()).toInt(), 1);

    // propanol -> methanol
    QCOMPARE(qRound(descriptor.value(methanol.get()).toDouble() * 100), 50);

    // propanol -> ethanol
    QCOMPARE(qRound(descriptor.value(ethanol.get()).toDouble() * 100), 75);
}

QTEST_APPLESS_MAIN(StructureSimilarityDescriptorTest)
