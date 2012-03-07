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

#include "pubchemwebtest.h"

#include <chemkit/pubchem.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

void PubChemWebTest::downloadFile()
{
    chemkit::PubChem pubchem;

    // CID 5950 is alanine
    chemkit::MoleculeFile *file = pubchem.downloadFile("5950");
    QVERIFY(file != 0);

    QCOMPARE(file->moleculeCount(), size_t(1));
    boost::shared_ptr<chemkit::Molecule> molecule = file->molecule();
    QCOMPARE(molecule->formula(), std::string("C3H7NO2"));

    delete file;
}

void PubChemWebTest::downloadMultiFile()
{
    chemkit::PubChem pubchem;

    QStringList ids;
    ids << "1" << "4" << "92" << "8" << "109" << "12";

    chemkit::MoleculeFile *file = pubchem.downloadFile(ids);
    QVERIFY(file != 0);

    QCOMPARE(file->moleculeCount(), size_t(6));

    for(int i = 0; i < ids.size(); i++){
        boost::shared_ptr<chemkit::Molecule> molecule = file->molecule(i);
        QVERIFY(molecule != 0);
        QByteArray id = ids[i].toAscii();
        QCOMPARE(molecule->name().c_str(), id.constData());
    }

    delete file;
}

void PubChemWebTest::search()
{
    chemkit::PubChem pubchem;

    // search for caffeine from is CAS number
    QStringList results = pubchem.search("58-08-2");
    QCOMPARE(results.size(), 1);
    QCOMPARE(results[0], QString("2519"));
}

void PubChemWebTest::standardizeFormula()
{
    chemkit::PubChem pubchem;

    std::string formula = pubchem.standardizeFormula("c3cccOc3", "smiles");
    QCOMPARE(formula, std::string("C1C=CC=CO1"));
}

QTEST_MAIN(PubChemWebTest)
