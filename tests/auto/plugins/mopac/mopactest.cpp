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

#include "mopactest.h"

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>
#include <chemkit/chemicalfileformat.h>

const std::string dataPath = "../../../data/";

void MopacTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::ChemicalFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "mopin") != formats.end());
    QVERIFY(std::find(formats.begin(), formats.end(), "mopcrt") != formats.end());
}

void MopacTest::ethanol()
{
    chemkit::ChemicalFile file(dataPath + "ethanol.mopin");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.moleculeCount(), 1);

    const chemkit::Molecule *molecule = file.molecule(0);
    QCOMPARE(molecule->formula(), std::string("C2H6O"));
}

void MopacTest::guanine()
{
    chemkit::ChemicalFile file(dataPath + "guanine.mopcrt");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.moleculeCount(), 1);

    const chemkit::Molecule *molecule = file.molecule(0);
    QCOMPARE(molecule->formula(), std::string("C5H5N5O"));
}

QTEST_APPLESS_MAIN(MopacTest)
