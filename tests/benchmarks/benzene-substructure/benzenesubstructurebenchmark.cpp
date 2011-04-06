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

// This benchmark measures the performance of the substructure
// isomorphism algorithms in chemkit.
//
// Based on: http://depth-first.com/articles/2009/01/22/mx-performance-comparison-3-substructure-search-in-mx-and-cdk

#include "benzenesubstructurebenchmark.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const std::string dataPath = "../../data/";

void BenzeneSubstructureBenchmark::benchmark()
{
    // load test file
    chemkit::MoleculeFile file(dataPath + "pubchem_416_benzenes.sdf");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // create benzene molecule
    chemkit::Molecule benzene("1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");

    QBENCHMARK {
        // number of substructure matches
        int count = 0;

        foreach(const chemkit::Molecule *molecule, file.molecules()){
            if(benzene.isSubstructureOf(molecule)){
                count++;
            }
        }

        // TODO: verify this number
        QCOMPARE(count, 412);
    }
}

QTEST_APPLESS_MAIN(BenzeneSubstructureBenchmark)
