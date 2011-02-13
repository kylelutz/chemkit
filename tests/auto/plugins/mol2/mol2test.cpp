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

class Mol2Test : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void read_data();
        void read();
};

void Mol2Test::initTestCase()
{
    QVERIFY(chemkit::ChemicalFileFormat::formats().contains("mol2"));
}

void Mol2Test::read_data()
{
    QTest::addColumn<QString>("fileName");
    QTest::addColumn<QString>("formula");

    QTest::newRow("uridine") << "uridine.mol2" << "C9H13N2O9P";
}

void Mol2Test::read()
{
    QFETCH(QString, fileName);
    QFETCH(QString, formula);

    chemkit::ChemicalFile file(dataPath + fileName);
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), formula);
}

QTEST_APPLESS_MAIN(Mol2Test)
#include "mol2test.moc"
