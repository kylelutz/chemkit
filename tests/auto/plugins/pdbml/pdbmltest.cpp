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

#include "pdbmltest.h"

#include <algorithm>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/polymerfileformat.h>

const std::string dataPath = "../../../data/";

void PdbmlTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::io::PolymerFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "pdbml") != formats.end());
}

void PdbmlTest::read_1UBQ()
{
    chemkit::io::PolymerFile file(dataPath + "1UBQ.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.polymerCount(), 1);

    // protein
    chemkit::Polymer *protein = file.polymer();
    QCOMPARE(protein->chainCount(), 1);

    // chain
    chemkit::PolymerChain *chain = protein->chain(0);
    QVERIFY(chain != 0);
    QCOMPARE(chain->residueCount(), 76);
    QCOMPARE(chain->sequenceString(), std::string("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                                  "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                                  "KESTLHLVLRLRGG"));

    // molecule
    QCOMPARE(protein->atomCount(), 660);
}

void PdbmlTest::read_2DHB()
{
    chemkit::io::PolymerFile file(dataPath + "2DHB.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.polymerCount(), 1);

    // protein
    chemkit::Polymer *protein = file.polymer();
    QCOMPARE(protein->chainCount(), 2);

    // chain A
    chemkit::PolymerChain *chainA = protein->chain(0);
    QCOMPARE(chainA->residueCount(), 141);
    QCOMPARE(chainA->sequenceString(), std::string("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                                   "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                                   "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                                   "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                                   "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::PolymerChain *chainB = protein->chain(1);
    QCOMPARE(chainB->residueCount(), 146);
    QCOMPARE(chainB->sequenceString(), std::string("VQLSGEEKAAVLALWDKVNEEEVGGEALGR"
                                                   "LLVVYPWTQRFFDSFGDLSNPGAVMGNPKV"
                                                   "KAHGKKVLHSFGEGVHHLDNLKGTFAALSE"
                                                   "LHCDKLHVDPENFRLLGNVLALVVARHFGK"
                                                   "DFTPELQASYQKVVAGVANALAHKYH"));
}

QTEST_APPLESS_MAIN(PdbmlTest)
