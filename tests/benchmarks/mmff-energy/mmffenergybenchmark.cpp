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

const std::string dataPath = "../../data/";

class MmffEnergyBenchmark : public QObject
{
    Q_OBJECT

    private slots:
        void benchmark();
};

void MmffEnergyBenchmark::benchmark()
{
    // load test file
    chemkit::ChemicalFile file(dataPath + "MMFF94_hypervalent.mol2");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::ForceField *mmff = chemkit::ForceField::create("mmff");
    mmff->setup();
    delete mmff;

    // total energy of all 753 molecules
    double totalEnergy = 0;

    QBENCHMARK_ONCE {
        foreach(const chemkit::Molecule *molecule, file.molecules()){
            chemkit::ForceField *forceField = chemkit::ForceField::create("mmff");
            QVERIFY(forceField);

            forceField->addMolecule(molecule);
            forceField->setup();
            //QVERIFY(forceField->isSetup());

            totalEnergy += forceField->energy();
        }
    }

    // expected total energy = 5228.05954
    QCOMPARE(qRound(totalEnergy), 5228);
}

QTEST_APPLESS_MAIN(MmffEnergyBenchmark)
#include "mmffenergybenchmark.moc"
