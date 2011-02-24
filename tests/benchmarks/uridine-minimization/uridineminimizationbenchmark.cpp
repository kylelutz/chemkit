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

#include <QtTest>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>
#include <chemkit/chemicalfile.h>

const QString dataPath = "../../data/";

class UridineMinimizationBenchmark : public QObject
{
    Q_OBJECT

    private slots:
        void benchmark();
};

void UridineMinimizationBenchmark::benchmark()
{
    chemkit::Molecule *molecule = chemkit::ChemicalFile::quickRead(dataPath + "uridine.mol2");
    QVERIFY(molecule != 0);

    chemkit::ForceField *forceField = chemkit::ForceField::create("uff");
    QVERIFY(forceField != 0);

    forceField->addMolecule(molecule);
    bool ok = forceField->setup();
    QVERIFY(ok);

    QBENCHMARK {
        for(;;){
            // converge when rmsg = 0.1
            bool converged = forceField->minimizationStep(0.1);

            if(converged){
                break;
            }
        }
    }

    delete forceField;
    delete molecule;
}

QTEST_APPLESS_MAIN(UridineMinimizationBenchmark)
#include "uridineminimizationbenchmark.moc"
