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

#include "wienerindextest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void WienerIndexTest::initTestCase()
{
    std::vector<std::string> descriptors = chemkit::MolecularDescriptor::descriptors();
    QVERIFY(std::find(descriptors.begin(), descriptors.end(), "wiener-index") != descriptors.end());
}

void WienerIndexTest::ethane()
{
    chemkit::Molecule ethane("CC", "smiles");
    QCOMPARE(ethane.formula(), std::string("C2H6"));
    QCOMPARE(ethane.descriptor("wiener-index").toInt(), 58);
}

QTEST_APPLESS_MAIN(WienerIndexTest)
