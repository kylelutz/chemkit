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

#include "rotatablebondstest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void RotatableBondsTest::initTestCase()
{
    std::vector<std::string> descriptors = chemkit::MolecularDescriptor::descriptors();
    QVERIFY(std::find(descriptors.begin(), descriptors.end(), "rotatable-bonds") != descriptors.end());
}

void RotatableBondsTest::test_data()
{
    QTest::addColumn<QString>("smiles");
    QTest::addColumn<int>("rotatableBonds");

    QTest::newRow("alanine") << "CC(C(=O)O)N" << 1;
    QTest::newRow("benzene") << "c1ccccc1" << 0;
    QTest::newRow("biphenyl") << "c1ccccc1(c2ccccc2)" << 1;
    QTest::newRow("isoleucine") << "CCC(C)C(C(=O)O)N" << 3;
    QTest::newRow("asparagine") << "C(C(C(=O)O)N)C(=O)N" << 3;
    QTest::newRow("octane") << "CCCCCCCC" << 5;
}

void RotatableBondsTest::test()
{
    QFETCH(QString, smiles);
    QFETCH(int, rotatableBonds);

    chemkit::Molecule molecule(smiles.toStdString(), "smiles");
    if(molecule.isEmpty()){
        qDebug() << "failed to load molecule";
    }
    QVERIFY(!molecule.isEmpty());

    QCOMPARE(molecule.descriptor("rotatable-bonds").toInt(), rotatableBonds);
}

QTEST_APPLESS_MAIN(RotatableBondsTest)
