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

#include "forcefieldtest.h" 

#include <algorithm>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>

#include "mockforcefield.h"

void ForceFieldTest::initTestCase()
{
    m_plugin = new MockForceFieldPlugin;

    std::vector<std::string> forceFields = chemkit::ForceField::forceFields();
    QVERIFY(std::find(forceFields.begin(), forceFields.end(), "mock") != forceFields.end());
}

void ForceFieldTest::create()
{
    chemkit::ForceField *forceField = chemkit::ForceField::create("mock");
    QVERIFY(forceField != 0);
    delete forceField;

    forceField = chemkit::ForceField::create("invalid_name");
    QVERIFY(forceField == 0);
}

void ForceFieldTest::name()
{
    chemkit::ForceField *forceField = chemkit::ForceField::create("mock");
    QCOMPARE(forceField->name(), std::string("mock"));
    delete forceField;
}

void ForceFieldTest::cleanupTestCase()
{
    delete m_plugin;

    std::vector<std::string> forceFields = chemkit::ForceField::forceFields();
    QVERIFY(std::find(forceFields.begin(), forceFields.end(), "mock") == forceFields.end());
}

QTEST_APPLESS_MAIN(ForceFieldTest)
