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
#include <chemkit/forcefield.h>
#include <chemkit/chemicalfile.h>

const QString dataPath = "../../../data/";

// The OplsTest class validates the OPLS force field implementation
// in the opls plugin. Energies were checked against those
// calculated by TINKER.
class OplsTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void energy_data();
        void energy();
};

void OplsTest::initTestCase()
{
    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "opls") != typers.end());

    QVERIFY(chemkit::ForceField::forceFields().contains("opls"));
}

void OplsTest::energy_data()
{
    QTest::addColumn<QString>("fileName");
    QTest::addColumn<QString>("formula");
    QTest::addColumn<double>("energy");

    QTest::newRow("water") << "water.mol" << "H2O" << 1.8698;
    QTest::newRow("methanol") << "methanol.sdf" << "CH4O" << 5.6693;
    QTest::newRow("ethanol") << "ethanol.cml" << "C2H6O" << 1.2309;

}

void OplsTest::energy()
{
    QFETCH(QString, fileName);
    QFETCH(QString, formula);
    QFETCH(double, energy);

    chemkit::Molecule *molecule = chemkit::ChemicalFile::quickRead(dataPath + fileName);
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->formula(), formula.toStdString());

    chemkit::ForceField *opls = chemkit::ForceField::create("opls");
    QVERIFY(opls != 0);

    opls->addMolecule(molecule);
    bool setup = opls->setup();
    QVERIFY(setup == true);

    QCOMPARE(qRound(opls->energy()), qRound(energy));

    delete molecule;
    delete opls;
}

QTEST_APPLESS_MAIN(OplsTest)
#include "oplstest.moc"
