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

#include <boost/range/algorithm.hpp>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/polymerfileformat.h>

const std::string dataPath = "../../../data/";

void PdbTest::initTestCase()
{
    // verify that the pdb plugin registered itself correctly
    QVERIFY(boost::count(chemkit::PolymerFileFormat::formats(), "pdb") == 1);
    QVERIFY(boost::count(chemkit::PolymerFileFormat::formats(), "pdbml") == 1);
}

void PdbTest::read_1BNA()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1BNA.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check nucleic acid
    QCOMPARE(file.polymerCount(), size_t(1));
    const boost::shared_ptr<chemkit::Polymer> &polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->name(),
             std::string("STRUCTURE OF A B-DNA DODECAMER. CONFORMATION AND DYNAMICS"));
    QCOMPARE(polymer->chainCount(), size_t(2));

    // check chains
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), size_t(12));
    QCOMPARE(chainA->sequenceString(), std::string("CGCGAATTCGCG"));

    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), size_t(12));
    QCOMPARE(chainB->sequenceString(), std::string("CGCGAATTCGCG"));

    // check ligands
    QCOMPARE(file.ligandCount(), size_t(80));
    foreach(const boost::shared_ptr<chemkit::Molecule> &ligand, file.ligands()){
        QCOMPARE(ligand->name(), std::string("HOH"));
        QCOMPARE(ligand->atomCount(), size_t(1));
        QCOMPARE(ligand->formula(), std::string("O"));
    }
}

void PdbTest::read_1UBQ()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1UBQ.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), size_t(1));
    const boost::shared_ptr<chemkit::Polymer> &polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->name(),
             std::string("STRUCTURE OF UBIQUITIN REFINED AT 1.8 ANGSTROMS RESOLUTION"));
    QCOMPARE(polymer->chainCount(), size_t(1));

    // check chain
    chemkit::PolymerChain *chain = polymer->chain(0);
    QVERIFY(chain != 0);

    // check residues
    QCOMPARE(chain->residueCount(), size_t(76));

    // check sequence
    QCOMPARE(chain->sequenceString(), std::string("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                                  "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                                  "KESTLHLVLRLRGG"));
}

void PdbTest::read_1UBQ_pdbml()
{
    chemkit::PolymerFile file(dataPath + "1UBQ.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.polymerCount(), size_t(1));

    // protein
    const boost::shared_ptr<chemkit::Polymer> &protein = file.polymer();
    QCOMPARE(protein->chainCount(), size_t(1));

    // chain
    chemkit::PolymerChain *chain = protein->chain(0);
    QVERIFY(chain != 0);
    QCOMPARE(chain->residueCount(), size_t(76));
    QCOMPARE(chain->sequenceString(), std::string("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                                  "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                                  "KESTLHLVLRLRGG"));

    // molecule
    QCOMPARE(protein->atomCount(), size_t(660));
}

void PdbTest::read_2DHB()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "2DHB.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), size_t(1));
    const boost::shared_ptr<chemkit::Polymer> &polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->name(),
             std::string("THREE DIMENSIONAL FOURIER SYNTHESIS OF HORSE "
                         "DEOXYHAEMOGLOBIN AT 2.8 ANGSTROMS RESOLUTION"));
    QCOMPARE(polymer->chainCount(), size_t(2));

    // chain A
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), size_t(141));
    QCOMPARE(chainA->sequenceString(), std::string("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                                   "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                                   "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                                   "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                                   "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), size_t(146));
    QCOMPARE(chainB->sequenceString(), std::string("VQLSGEEKAAVLALWDKVNEEEVGGEALGR"
                                                   "LLVVYPWTQRFFDSFGDLSNPGAVMGNPKV"
                                                   "KAHGKKVLHSFGEGVHHLDNLKGTFAALSE"
                                                   "LHCDKLHVDPENFRLLGNVLALVVARHFGK"
                                                   "DFTPELQASYQKVVAGVANALAHKYH"));

    // check ligands
    QCOMPARE(file.ligandCount(), size_t(4));
    QCOMPARE(file.ligand(0)->name(), std::string("PROTOPORPHYRIN IX CONTAINING FE"));
    QCOMPARE(file.ligand(1)->name(), std::string("PROTOPORPHYRIN IX CONTAINING FE"));
    QCOMPARE(file.ligand(2)->name(), std::string("HOH"));
    QCOMPARE(file.ligand(3)->name(), std::string("HOH"));
}

void PdbTest::read_2DHB_pdbml()
{
    chemkit::PolymerFile file(dataPath + "2DHB.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(file.polymerCount(), size_t(1));

    // protein
    const boost::shared_ptr<chemkit::Polymer> &protein = file.polymer();
    QCOMPARE(protein->chainCount(), size_t(2));

    // chain A
    chemkit::PolymerChain *chainA = protein->chain(0);
    QCOMPARE(chainA->residueCount(), size_t(141));
    QCOMPARE(chainA->sequenceString(), std::string("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                                   "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                                   "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                                   "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                                   "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::PolymerChain *chainB = protein->chain(1);
    QCOMPARE(chainB->residueCount(), size_t(146));
    QCOMPARE(chainB->sequenceString(), std::string("VQLSGEEKAAVLALWDKVNEEEVGGEALGR"
                                                   "LLVVYPWTQRFFDSFGDLSNPGAVMGNPKV"
                                                   "KAHGKKVLHSFGEGVHHLDNLKGTFAALSE"
                                                   "LHCDKLHVDPENFRLLGNVLALVVARHFGK"
                                                   "DFTPELQASYQKVVAGVANALAHKYH"));
}

void PdbTest::read_alphabet()
{
    // create file
    chemkit::PolymerFile file(dataPath + "alphabet.pdb");

    // read file
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), size_t(1));
    const boost::shared_ptr<chemkit::Polymer> &polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), size_t(1));

    chemkit::PolymerChain *chain = polymer->chain(0);
    QCOMPARE(chain->residueCount(), size_t(20));
    QCOMPARE(chain->sequenceString(), std::string("ADNRCEQGHILKMFPSTWYV"));
}

void PdbTest::read_fmc()
{
    chemkit::PolymerFile file(dataPath + "fmc.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.polymerCount(), size_t(0));
    QCOMPARE(file.ligandCount(), size_t(1));

    const boost::shared_ptr<chemkit::Molecule> &molecule = file.ligand(0);
    QCOMPARE(molecule->atomCount(), size_t(2596));
}

QTEST_APPLESS_MAIN(PdbTest)
