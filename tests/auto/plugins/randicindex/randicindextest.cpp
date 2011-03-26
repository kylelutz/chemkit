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
#include <chemkit/moleculardescriptor.h>

class RandicIndexTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void ethane();
        void isobutane();
        void dimethylpropane();
        void octane();
};

void RandicIndexTest::initTestCase()
{
    QVERIFY(chemkit::MolecularDescriptor::descriptors().contains("randic-index"));
}

void RandicIndexTest::ethane()
{
    chemkit::Molecule ethane("CC", "smiles");
    QCOMPARE(ethane.formula(), std::string("C2H6"));

    // index = 1.0
    QCOMPARE(ethane.descriptor("randic-index").toInt(), 1);
}

void RandicIndexTest::isobutane()
{
    chemkit::Molecule isobutane("CC(C)C", "smiles");
    QCOMPARE(isobutane.formula(), std::string("C4H10"));

    // index = 1.7321
    QCOMPARE(qRound(isobutane.descriptor("randic-index").toDouble()), 2);
}

void RandicIndexTest::dimethylpropane()
{
    chemkit::Molecule dimethylpropane("CC(C)(C)C", "smiles");
    QCOMPARE(dimethylpropane.formula(), std::string("C5H12"));

    // index = 2.0
    QCOMPARE(qRound(dimethylpropane.descriptor("randic-index").toDouble()), 2);
}

void RandicIndexTest::octane()
{
    chemkit::Molecule octane("CCCCCCCC", "smiles");
    QCOMPARE(octane.formula(), std::string("C8H18"));

    // index = 3.9142
    QCOMPARE(qRound(octane.descriptor("randic-index").toDouble()), 4);
}

QTEST_APPLESS_MAIN(RandicIndexTest)
#include "randicindextest.moc"
