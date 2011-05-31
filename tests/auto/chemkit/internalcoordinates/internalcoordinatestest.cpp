/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "internalcoordinatestest.h"

#include <chemkit/chemkit.h>
#include <chemkit/internalcoordinates.h>
#include <vector>

void InternalCoordinatesTest::size()
{
    chemkit::InternalCoordinates coordinates(1);
    QCOMPARE(coordinates.size(),1);

    chemkit::InternalCoordinates newCoordinates(coordinates);
    QCOMPARE(newCoordinates.size(),1);
}

void InternalCoordinatesTest::coordinates()
{
    chemkit::InternalCoordinates newCoordinates(1);
    QCOMPARE(newCoordinates.size(),1);

    chemkit::Float r = 0.0;
    chemkit::Float theta = 1.0;
    chemkit::Float phi = 2.0;
    newCoordinates.setCoordinates(0, r, theta, phi);
    std::vector<chemkit::Float> coordinates = newCoordinates.coordinates(0);

    QCOMPARE(coordinates[0], chemkit::Float(0.0));
    QCOMPARE(coordinates[1], chemkit::Float(1.0));
    QCOMPARE(coordinates[2], chemkit::Float(2.0));
}

void InternalCoordinatesTest::coordinatesRadians()
{
    chemkit::InternalCoordinates newCoordinates(1);
    QCOMPARE(newCoordinates.size(),1);

    chemkit::Float r = 0.0;
    chemkit::Float theta = 1.0;
    chemkit::Float phi = 2.0;
    newCoordinates.setCoordinatesRadians(0, r, theta, phi);

    r = r * chemkit::constants::RadiansToDegrees;
    theta = theta * chemkit::constants::RadiansToDegrees;
    phi = phi * chemkit::constants::RadiansToDegrees;
    std::vector<chemkit::Float> coordinates = newCoordinates.coordinates(0);
    QCOMPARE(coordinates[0], r);
    QCOMPARE(coordinates[1], theta);
    QCOMPARE(coordinates[2], phi);
}

void InternalCoordinatesTest::connections()
{
    chemkit::InternalCoordinates newCoordinates(1);
    QCOMPARE(newCoordinates.size(),1);

    newCoordinates.setConnections(0, 1, 2, 3);
    std::vector<int> connections = newCoordinates.connections(0);

    QCOMPARE(connections[0], 1);
    QCOMPARE(connections[1], 2);
    QCOMPARE(connections[2], 3);
}

QTEST_APPLESS_MAIN(InternalCoordinatesTest)
