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

#include "pdbtest.h"

#include <algorithm>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/polymerfileformat.h>

const std::string dataPath = "../../../data/";

void PdbTest::initTestCase()
{
    std::vector<std::string> formats = chemkit::io::PolymerFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "pdb") != formats.end());
}

void PdbTest::read_1BNA()
{
    // create file
    chemkit::io::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1BNA.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check nucleic acid
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 2);

    // check chains
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), 12);
    QCOMPARE(chainA->sequenceString(), std::string("CGCGAATTCGCG"));

    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), 12);
    QCOMPARE(chainB->sequenceString(), std::string("CGCGAATTCGCG"));
}

void PdbTest::read_1UBQ()
{
    // create file
    chemkit::io::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1UBQ.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 1);

    // check chain
    chemkit::PolymerChain *chain = polymer->chain(0);
    QVERIFY(chain != 0);

    // check residues
    QCOMPARE(chain->residueCount(), 76);

    // check sequence
    QCOMPARE(chain->sequenceString(), std::string("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                                  "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                                  "KESTLHLVLRLRGG"));
}

void PdbTest::read_2DHB()
{
    // create file
    chemkit::io::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "2DHB.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 2);

    // chain A
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), 141);
    QCOMPARE(chainA->sequenceString(), std::string("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                                   "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                                   "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                                   "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                                   "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), 146);
    QCOMPARE(chainB->sequenceString(), std::string("VQLSGEEKAAVLALWDKVNEEEVGGEALGR"
                                                   "LLVVYPWTQRFFDSFGDLSNPGAVMGNPKV"
                                                   "KAHGKKVLHSFGEGVHHLDNLKGTFAALSE"
                                                   "LHCDKLHVDPENFRLLGNVLALVVARHFGK"
                                                   "DFTPELQASYQKVVAGVANALAHKYH"));
}

void PdbTest::read_alphabet()
{
    // create file
    chemkit::io::PolymerFile file(dataPath + "alphabet.pdb");

    // read file
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 1);

    chemkit::PolymerChain *chain = polymer->chain(0);
    QCOMPARE(chain->residueCount(), 20);
    QCOMPARE(chain->sequenceString(), std::string("ADNRCEQGHILKMFPSTWYV"));
}

void PdbTest::read_fmc()
{
    chemkit::io::PolymerFile file(dataPath + "fmc.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);
}

QTEST_APPLESS_MAIN(PdbTest)
