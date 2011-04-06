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

#include "ufftest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>

void UffTest::initTestCase()
{
    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "uff") != typers.end());

    std::vector<std::string> forceFields = chemkit::ForceField::forceFields();
    QVERIFY(std::find(forceFields.begin(), forceFields.end(), "uff") != forceFields.end());
}

QTEST_APPLESS_MAIN(UffTest)
