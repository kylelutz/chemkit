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

// This benchmark measures the time it takes to calculate the
// solvent accessible surface area of the protein hemoglobin
// (PDB ID: 2DHB). The protein contains 146 residues and 2201
// atoms.

#include <QtTest>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/molecularsurface.h>

const std::string dataPath = "../../data/";

class ProteinSurfaceBenchmark : public QObject
{
    Q_OBJECT

    private slots:
        void benchmark();
};

void ProteinSurfaceBenchmark::benchmark()
{
    chemkit::PolymerFile file(dataPath + "2DHB.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 2201);

    QBENCHMARK {
        chemkit::MolecularSurface surface(protein);
        surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);

        QCOMPARE(qRound(surface.surfaceArea()), 14791);
    }
}

QTEST_APPLESS_MAIN(ProteinSurfaceBenchmark)
#include "proteinsurfacebenchmark.moc"
