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

#include "polymertest.h" 

#include <chemkit/polymer.h>

void PolymerTest::chains()
{
    chemkit::Polymer polymer;
    chemkit::PolymerChain *chain1 = polymer.addChain();
    QCOMPARE(polymer.chainCount(), 1);
    QVERIFY(polymer.chain(0) == chain1);
    QVERIFY(polymer.chains()[0] == chain1);

    chemkit::PolymerChain *chain2 = polymer.addChain();
    QCOMPARE(polymer.chainCount(), 2);
    QVERIFY(polymer.chain(1) == chain2);
    QVERIFY(polymer.chains()[0] == chain1);
    QVERIFY(polymer.chains()[1] == chain2);

    polymer.removeChain(chain1);
    QCOMPARE(polymer.chainCount(), 1);
    QVERIFY(polymer.chain(0) == chain2);
}

QTEST_APPLESS_MAIN(PolymerTest)
