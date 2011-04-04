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

#include "chemicalfiletest.h"

#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>

void ChemicalFileTest::fileName()
{
    chemkit::ChemicalFile file;
    QCOMPARE(file.fileName(), std::string());

    file.setFileName("foo");
    QCOMPARE(file.fileName(), std::string("foo"));

    file.setFileName("bar");
    QCOMPARE(file.fileName(), std::string("bar"));

    file.setFileName(std::string());
    QCOMPARE(file.fileName(), std::string());

    chemkit::ChemicalFile fileWithName("foobar");
    QCOMPARE(fileWithName.fileName(), std::string("foobar"));
}

void ChemicalFileTest::format()
{
    chemkit::ChemicalFile file;
    QVERIFY(file.format() == 0);
}

void ChemicalFileTest::contains()
{
    chemkit::ChemicalFile file;
    QCOMPARE(file.contains(0), false);

    chemkit::Molecule *molecule = new chemkit::Molecule;
    file.addMolecule(molecule);
    QCOMPARE(file.contains(molecule), true);

    chemkit::Molecule *anotherMolecule = new chemkit::Molecule;
    QCOMPARE(file.contains(anotherMolecule), false);

    file.addMolecule(anotherMolecule);
    QCOMPARE(file.contains(anotherMolecule), true);

    file.removeMolecule(molecule);
    QCOMPARE(file.contains(molecule), false);

    file.addMolecule(molecule);
}

void ChemicalFileTest::fileData()
{
    chemkit::ChemicalFile file;

    file.setFileData("foo", "bar");
    QCOMPARE(file.fileData("foo").toString(), QString("bar"));

    file.setFileData("number", 4);
    QCOMPARE(file.fileData("number").toInt(), 4);
}

void ChemicalFileTest::moleculeData()
{
    chemkit::ChemicalFile file;
    chemkit::Molecule *molecule = new chemkit::Molecule;
    file.addMolecule(molecule);

    file.setMoleculeData(molecule, "boilingPoint", 38);
    QCOMPARE(file.moleculeData(molecule, "boilingPoint").toInt(), 38);

    chemkit::Molecule *otherMolecule = new chemkit::Molecule;
    file.addMolecule(otherMolecule);
    file.setMoleculeData(otherMolecule, "boilingPoint", 87);
    QCOMPARE(file.moleculeData(otherMolecule, "boilingPoint").toInt(), 87);
    QCOMPARE(file.moleculeData(molecule, "boilingPoint").toInt(), 38);
}

QTEST_APPLESS_MAIN(ChemicalFileTest)
