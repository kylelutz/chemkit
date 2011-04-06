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

// Volume and surface area measurements for van der waals and solvent
// accessible surfaces were validated against the output of the 'asv'
// program. (http://petitjeanmichel.free.fr/itoweb.petitjean.freeware.html#ASV)

#include "molecularsurfacetest.h"

#include <chemkit/point3.h>
#include <chemkit/polymer.h>
#include <chemkit/vector3.h>
#include <chemkit/molecule.h>
#include <chemkit/polymerfile.h>
#include <chemkit/moleculefile.h>
#include <chemkit/molecularsurface.h>

const std::string dataPath = "../../../data/";

void MolecularSurfaceTest::molecule()
{
    chemkit::Molecule molecule;
    chemkit::MolecularSurface surface(&molecule);
    QVERIFY(surface.molecule() == &molecule);

    chemkit::MolecularSurface emptySurface;
    QVERIFY(emptySurface.molecule() == 0);
}

void MolecularSurfaceTest::probeRadius()
{
    chemkit::Molecule molecule;
    chemkit::MolecularSurface surface(&molecule);

    // ensure default probe radius is 1.4
    QCOMPARE(surface.probeRadius(), chemkit::Float(1.4));

    surface.setProbeRadius(2.5);
    QCOMPARE(surface.probeRadius(), chemkit::Float(2.5));

    surface.setProbeRadius(0);
    QCOMPARE(surface.probeRadius(), chemkit::Float(0.0));
}

void MolecularSurfaceTest::surfaceType()
{
    chemkit::Molecule molecule;
    chemkit::MolecularSurface surface(&molecule);

    // ensure default surface type is van der waals
    QCOMPARE(surface.surfaceType(), chemkit::MolecularSurface::VanDerWaals);

    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(surface.surfaceType(), chemkit::MolecularSurface::SolventAccessible);

    surface.setSurfaceType(chemkit::MolecularSurface::SolventExcluded);
    QCOMPARE(surface.surfaceType(), chemkit::MolecularSurface::SolventExcluded);

    // check surface type when set with the constructor
    chemkit::MolecularSurface surface2(&molecule, chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(surface2.surfaceType(), chemkit::MolecularSurface::SolventAccessible);
}

void MolecularSurfaceTest::hydrogen()
{
    chemkit::Molecule molecule;
    molecule.addAtom("H");

    chemkit::MolecularSurface surface(&molecule);
    QCOMPARE(qRound(surface.volume()), 7);
    QCOMPARE(qRound(surface.surfaceArea()), 18);

    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 74);
    QCOMPARE(qRound(surface.surfaceArea()), 85);

    chemkit::Atom *H2 = molecule.addAtom("H");
    H2->setPosition(2.4, 0, 0);
    surface.setMolecule(&molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 14);
    QCOMPARE(qRound(surface.surfaceArea()), 36);
}

void MolecularSurfaceTest::water()
{
}

void MolecularSurfaceTest::serine()
{
    chemkit::MoleculeFile file(dataPath + "serine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule);
    QCOMPARE(molecule->size(), 14);

    // van der waals surface
    chemkit::MolecularSurface surface(molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 93);
    QCOMPARE(qRound(surface.surfaceArea()), 129);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 360);
    QCOMPARE(qRound(surface.surfaceArea()), 262);
}

void MolecularSurfaceTest::guanine()
{
    chemkit::MoleculeFile file(dataPath + "guanine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule);
    QCOMPARE(molecule->size(), 16);

    // van der waals surface
    chemkit::MolecularSurface surface(molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 119);
    QCOMPARE(qRound(surface.surfaceArea()), 153);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 437);
    QCOMPARE(qRound(surface.surfaceArea()), 309);
}

void MolecularSurfaceTest::methane()
{
    chemkit::MoleculeFile file(dataPath + "methane.xyz");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule);
    QCOMPARE(molecule->size(), 5);

    // van der waals surface
    chemkit::MolecularSurface surface(molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 28);
    QCOMPARE(qRound(surface.surfaceArea()), 48);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 157);
    QCOMPARE(qRound(surface.surfaceArea()), 144);
}

void MolecularSurfaceTest::ethanol()
{
    chemkit::MoleculeFile file(dataPath + "ethanol.cml");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule);
    QCOMPARE(molecule->size(), 9);

    // van der waals surface
    chemkit::MolecularSurface surface(molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 53);
    QCOMPARE(qRound(surface.surfaceArea()), 81);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 245);
    QCOMPARE(qRound(surface.surfaceArea()), 199);
}

void MolecularSurfaceTest::adenosine()
{
    chemkit::MoleculeFile file(dataPath + "adenosine.mol");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule);
    QCOMPARE(molecule->size(), 32);

    // van der waals surface
    chemkit::MolecularSurface surface(molecule);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 219);
    QCOMPARE(qRound(surface.surfaceArea()), 273);

    // solvent accessible surface (probe radius = 1.4)
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    surface.setProbeRadius(1.4);
    QCOMPARE(qRound(surface.volume()), 723);
    QCOMPARE(qRound(surface.surfaceArea()), 457);

    // solvent accessible surface (probe radius = 1.0)
    surface.setProbeRadius(1.0);
    QCOMPARE(qRound(surface.volume()), 552);
    QCOMPARE(qRound(surface.surfaceArea()), 397);
}

void MolecularSurfaceTest::lysozyme()
{
    chemkit::PolymerFile file(dataPath + "2LYZ.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 1001);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 10729);
    QCOMPARE(qRound(surface.surfaceArea()), 12660);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 23737);
    QCOMPARE(qRound(surface.surfaceArea()), 6696);
}

void MolecularSurfaceTest::cytochrome()
{
    chemkit::PolymerFile file(dataPath + "3CYT.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 1600);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 17285);
    QCOMPARE(qRound(surface.surfaceArea()), 20616);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 39079);
    QCOMPARE(qRound(surface.surfaceArea()), 11533);
}

void MolecularSurfaceTest::toxin()
{
    chemkit::PolymerFile file(dataPath + "2SN3.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 948);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 6187);
    QCOMPARE(qRound(surface.surfaceArea()), 7214);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 13396);
    QCOMPARE(qRound(surface.surfaceArea()), 4383);
}

void MolecularSurfaceTest::hydrolase()
{
    chemkit::PolymerFile file(dataPath + "1THM.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 2003);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 21517);
    QCOMPARE(qRound(surface.surfaceArea()), 25741);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 44795);
    QCOMPARE(qRound(surface.surfaceArea()), 9852);
}

void MolecularSurfaceTest::hemoglobin()
{
    chemkit::PolymerFile file(dataPath + "2DHB.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 2201);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 24050);
    QCOMPARE(qRound(surface.surfaceArea()), 28917);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 54214);
    QCOMPARE(qRound(surface.surfaceArea()), 14791);
}

void MolecularSurfaceTest::dna()
{
    chemkit::PolymerFile file(dataPath + "1BNA.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *nucleicAcid = file.polymer();
    QVERIFY(nucleicAcid);
    QCOMPARE(nucleicAcid->size(), 486);

    // van der waals surface
    chemkit::MolecularSurface surface(nucleicAcid);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 4845);
    QCOMPARE(qRound(surface.surfaceArea()), 5615);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 11545);
    QCOMPARE(qRound(surface.surfaceArea()), 4660);
}

void MolecularSurfaceTest::ribozyme()
{
    chemkit::PolymerFile file(dataPath + "1MME.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 1746);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 17103);
    QCOMPARE(qRound(surface.surfaceArea()), 19404);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 38918);
    QCOMPARE(qRound(surface.surfaceArea()), 13523);
}

void MolecularSurfaceTest::ubiqutin()
{
    chemkit::PolymerFile file(dataPath + "1UBQ.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    chemkit::Polymer *protein = file.polymer();
    QVERIFY(protein);
    QCOMPARE(protein->size(), 602);

    // van der waals surface
    chemkit::MolecularSurface surface(protein);
    surface.setSurfaceType(chemkit::MolecularSurface::VanDerWaals);
    QCOMPARE(qRound(surface.volume()), 6559);
    QCOMPARE(qRound(surface.surfaceArea()), 7916);

    // solvent accessible surface
    surface.setSurfaceType(chemkit::MolecularSurface::SolventAccessible);
    QCOMPARE(qRound(surface.volume()), 15414);
    QCOMPARE(qRound(surface.surfaceArea()), 4871);
}

QTEST_APPLESS_MAIN(MolecularSurfaceTest)
