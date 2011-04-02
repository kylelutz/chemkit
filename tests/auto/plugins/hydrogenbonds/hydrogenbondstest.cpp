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
#include <chemkit/moleculardescriptor.h>

class HydrogenBondsTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void ethanol();
        void guanine();
};

void HydrogenBondsTest::initTestCase()
{
    std::vector<std::string> descriptors = chemkit::MolecularDescriptor::descriptors();
    QVERIFY(std::find(descriptors.begin(), descriptors.end(), "hydrogen-bond-donors") != descriptors.end());
    QVERIFY(std::find(descriptors.begin(), descriptors.end(), "hydrogen-bond-acceptors") != descriptors.end());
}

void HydrogenBondsTest::ethanol()
{
    chemkit::Molecule ethanol("CCO", "smiles");
    QCOMPARE(ethanol.formula(), std::string("C2H6O"));

    QCOMPARE(ethanol.descriptor("hydrogen-bond-donors").toInt(), 1);
    QCOMPARE(ethanol.descriptor("hydrogen-bond-acceptors").toInt(), 1);
}

void HydrogenBondsTest::guanine()
{
    chemkit::Molecule guanine("c1[nH]c2c(n1)c(=O)[nH]c(n2)N", "smiles");
    QCOMPARE(guanine.formula(), std::string("C5H5N5O"));

    QCOMPARE(guanine.descriptor("hydrogen-bond-donors").toInt(), 4);
    QCOMPARE(guanine.descriptor("hydrogen-bond-acceptors").toInt(), 6);
}

QTEST_APPLESS_MAIN(HydrogenBondsTest)
#include "hydrogenbondstest.moc"
