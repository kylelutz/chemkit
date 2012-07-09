/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "coordinatepredictortest.h"

#include <limits>

#include <chemkit/molecule.h>
#include <chemkit/coordinatepredictor.h>

void CoordinatePredictorTest::molecule()
{
    chemkit::CoordinatePredictor predictor;
    QVERIFY(predictor.molecule() == 0);

    chemkit::Molecule molecule;
    predictor.setMolecule(&molecule);
    QVERIFY(predictor.molecule() == &molecule);

    predictor.setMolecule(0);
    QVERIFY(predictor.molecule() == 0);
}

void CoordinatePredictorTest::eliminateCloseContacts()
{
    // create ethanol molecule
    chemkit::Molecule ethanol("CCO", "smiles");
    QCOMPARE(ethanol.formula(), std::string("C2H6O"));
    QVERIFY(ethanol.distance(ethanol.atom(0), ethanol.atom(1)) == chemkit::Real(0));

    // eliminate all close atom contacts less than two angstroms
    bool modified = chemkit::CoordinatePredictor::eliminateCloseContacts(&ethanol, 2.0);
    QVERIFY(modified == true);
    QVERIFY(ethanol.distance(ethanol.atom(0), ethanol.atom(1)) != chemkit::Real(0));

    chemkit::Real closestDistance = std::numeric_limits<chemkit::Real>::max();
    for(size_t i = 0; i < ethanol.size(); i++){
        for(size_t j = i + 1; j < ethanol.size(); j++){
            chemkit::Real distance = ethanol.distance(ethanol.atom(i),
                                                      ethanol.atom(j));

            if(distance < closestDistance){
                closestDistance = distance;
            }
        }
    }

    // verify that no two atoms are less than two angstroms from each other
    QVERIFY(closestDistance >= 2.0);

    // run the algorithm again but verify that nothing is modifed
    modified = chemkit::CoordinatePredictor::eliminateCloseContacts(&ethanol);
    QVERIFY(modified == false);
}

QTEST_APPLESS_MAIN(CoordinatePredictorTest)
