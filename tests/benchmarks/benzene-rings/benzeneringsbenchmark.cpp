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

// This benchmark measures the performance of the ring perception
// algorithm.
//
// Based on: http://depth-first.com/articles/2009/01/21/mx-performance-comparison-2-exhaustive-ring-perception-in-mx-and-cdk

#include "benzeneringsbenchmark.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const std::string dataPath = "../../data/";

void BenzeneRingsBenchmark::benchmark()
{
    // load test file
    chemkit::MoleculeFile file(dataPath + "pubchem_416_benzenes.sdf");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // total number of rings in all 416 molecules
    int ringCount = 0;

    QBENCHMARK_ONCE {
        foreach(const chemkit::Molecule *molecule, file.molecules()){
            // don't use molecule->ringCount() because it
            // may not actually perceive the rings.
            ringCount += molecule->rings().size();
        }
    }

    // the ringCount in the MX library is 2179 because they
    // calculate the exhaustive set of rings instead of the
    // SSSR like chemkit does.
    //
    // TODO: verify this number
    QCOMPARE(ringCount, 1288);
}

QTEST_APPLESS_MAIN(BenzeneRingsBenchmark)
