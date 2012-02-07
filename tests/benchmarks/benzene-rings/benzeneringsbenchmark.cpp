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

// This benchmark measures the performance of the ring perception
// algorithm.
//
// Based on: http://depth-first.com/articles/2009/01/21/mx-performance-comparison-2-exhaustive-ring-perception-in-mx-and-cdk

#include "benzeneringsbenchmark.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

const std::string dataPath = "../../data/";

void BenzeneRingsBenchmark::benchmark()
{
    // load test file
    chemkit::MoleculeFile file(dataPath + "pubchem_416_benzenes.sdf");
    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);

    // total number of rings in all 416 molecules
    int ringCount = 0;

    QBENCHMARK_ONCE {
        foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, file.molecules()){
            // don't use molecule->ringCount() because it
            // may not actually perceive the rings.
            ringCount += molecule->rings().size();
        }
    }

    // the ringCount in the MX library is 2179 because they
    // calculate the exhaustive set of rings instead of the
    // SSSR like chemkit does.
    //
    // TODO: verify this number
    QCOMPARE(ringCount, 1288);
}

QTEST_APPLESS_MAIN(BenzeneRingsBenchmark)
