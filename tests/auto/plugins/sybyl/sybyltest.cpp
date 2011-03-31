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

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/chemicalfile.h>
#include <chemkit/chemicalfileformat.h>

const std::string dataPath = "../../../data/";

class SybylTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void readMol2_data();
        void readMol2();
};

void SybylTest::initTestCase()
{
    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "sybyl") != typers.end());

    std::vector<std::string> fileFormats = chemkit::ChemicalFileFormat::formats();
    QVERIFY(std::find(fileFormats.begin(), fileFormats.end(), "mol2") != fileFormats.end());
}

void SybylTest::readMol2_data()
{
    QTest::addColumn<QString>("fileName");
    QTest::addColumn<QString>("formula");

    QTest::newRow("uridine") << "uridine.mol2" << "C9H13N2O9P";
}

void SybylTest::readMol2()
{
    QFETCH(QString, fileName);
    QFETCH(QString, formula);

    chemkit::ChemicalFile file(dataPath + fileName.toStdString());
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), 1);
    chemkit::Molecule *molecule = file.molecule();
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), formula.toStdString());
}

QTEST_APPLESS_MAIN(SybylTest)
#include "sybyltest.moc"
