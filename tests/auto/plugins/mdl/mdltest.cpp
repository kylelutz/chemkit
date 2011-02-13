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
#include <chemkit/chemicalfile.h>
#include <chemkit/chemicalfileformat.h>

const QString dataPath = "../../../data/";

class MdlTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();

        void read_methanol();
        void read_guanine();
        void read_benzenes();
};

void MdlTest::initTestCase()
{
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("mol"));
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("mdl"));
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("sdf"));
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("sd"));
}

void MdlTest::read_methanol()
{
    chemkit::ChemicalFile file(dataPath + "methanol.sdf");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    // check molecule
    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QCOMPARE(molecule->formula(), QString("CH4O"));

    // check data
    QCOMPARE(molecule->name(), QString("887"));
    QCOMPARE(file.moleculeData(molecule, "PUBCHEM_COMPOUND_CID").toString(), QString("887"));
    QCOMPARE(file.moleculeData(molecule, "PUBCHEM_HEAVY_ATOM_COUNT").toInt(), 2);
}

void MdlTest::read_guanine()
{
    chemkit::ChemicalFile file(dataPath + "guanine.mol");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QCOMPARE(ok, true);

    // check format
    QVERIFY(file.format() != 0);
    QCOMPARE(file.formatName(), QString("mol"));

    // check molecule
    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *guanine = file.molecule();
    QCOMPARE(guanine->formula(), QString("C5H5N5O"));
    QCOMPARE(guanine->name(), QString("Guanine"));
    QCOMPARE(guanine->atomCount(), 16);
    QCOMPARE(guanine->bondCount(), 17);
}

void MdlTest::read_benzenes()
{
    chemkit::ChemicalFile file(dataPath + "pubchem_416_benzenes.sdf");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QCOMPARE(ok, true);

    // check format
    QVERIFY(file.format() != 0);
    QCOMPARE(file.formatName(), QString("sdf"));

    // check molecules
    QCOMPARE(file.moleculeCount(), 416);

    // check molecule data
    foreach(const chemkit::Molecule *molecule, file.molecules()){
        QCOMPARE(molecule->name(), file.moleculeData(molecule, "PUBCHEM_COMPOUND_CID").toString());
    }
}

QTEST_APPLESS_MAIN(MdlTest)
#include "mdltest.moc"
