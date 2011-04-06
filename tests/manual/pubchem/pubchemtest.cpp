/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include "pubchemtest.h"

#include <chemkit/pubchem.h>
#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>

void PubChemTest::downloadFile()
{
    chemkit::PubChem pubchem;

    // CID 5950 is alanine
    chemkit::ChemicalFile *file = pubchem.downloadFile("5950");
    QVERIFY(file != 0);

    QCOMPARE(file->moleculeCount(), 1);
    chemkit::Molecule *molecule = file->molecule();
    QCOMPARE(molecule->formula(), std::string("C3H7NO2"));

    delete file;
}

void PubChemTest::downloadMultiFile()
{
    chemkit::PubChem pubchem;

    QStringList ids;
    ids << "1" << "4" << "92" << "8" << "109" << "12";

    chemkit::ChemicalFile *file = pubchem.downloadFile(ids);
    QVERIFY(file != 0);

    QCOMPARE(file->moleculeCount(), 6);

    for(int i = 0; i < ids.size(); i++){
        chemkit::Molecule *molecule = file->molecule(i);
        QVERIFY(molecule != 0);
        QCOMPARE(molecule->name(), ids[i].toStdString());
    }

    delete file;
}

void PubChemTest::search()
{
    chemkit::PubChem pubchem;

    // search for caffeine from is CAS number
    QStringList results = pubchem.search("58-08-2");
    QCOMPARE(results.size(), 1);
    QCOMPARE(results[0], QString("2519"));
}

void PubChemTest::standardizeFormula()
{
    chemkit::PubChem pubchem;

    std::string formula = pubchem.standardizeFormula("c3cccOc3", "smiles");
    QCOMPARE(formula, std::string("C1C=CC=CO1"));
}

QTEST_MAIN(PubChemTest)
