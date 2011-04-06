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

// This benchmark reads a 33 molecule sdf file and calculates
// the molecular masses for each molecule.
//
// Based on: http://depth-first.com/articles/2009/01/20/open-benchmarks-for-cheminformatics-first-performance-comparison-between-cdk-and-mx

#include "molecularmassesbenchmark.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

void MolecularMassesBenchmark::benchmark()
{
    QBENCHMARK {
        chemkit::MoleculeFile file("pubchem_sample_33.sdf");
        bool ok = file.read();
        if(!ok)
            qDebug() << file.errorString().c_str();
        QVERIFY(ok);

        double totalMass = 0;

        foreach(const chemkit::Molecule *molecule, file.molecules()){
            totalMass += molecule->mass();
        }

        QCOMPARE(qRound(totalMass), 6799);
    }
}

QTEST_APPLESS_MAIN(MolecularMassesBenchmark)
