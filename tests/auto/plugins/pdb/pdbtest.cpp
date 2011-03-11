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

#include <QtTest>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/polymerfileformat.h>

const QString dataPath = "../../../data/";

class PdbTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void read_1BNA();
        void read_1UBQ();
        void read_2DHB();
        void read_alphabet();
        void read_fmc();
};

void PdbTest::initTestCase()
{
    QVERIFY(chemkit::PolymerFileFormat::formats().contains("pdb"));
}

void PdbTest::read_1BNA()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1BNA.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check nucleic acid
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 2);

    // check chains
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), 12);
    QCOMPARE(chainA->sequenceString(), QString("CGCGAATTCGCG"));

    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), 12);
    QCOMPARE(chainB->sequenceString(), QString("CGCGAATTCGCG"));
}

void PdbTest::read_1UBQ()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "1UBQ.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
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
    QCOMPARE(chain->sequenceString(), QString("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                              "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                              "KESTLHLVLRLRGG"));
}

void PdbTest::read_2DHB()
{
    // create file
    chemkit::PolymerFile file;

    // read file
    bool ok = file.read(dataPath + "2DHB.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 2);

    // chain A
    chemkit::PolymerChain *chainA = polymer->chain(0);
    QCOMPARE(chainA->residueCount(), 141);
    QCOMPARE(chainA->sequenceString(), QString("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                               "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                               "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                               "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                               "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::PolymerChain *chainB = polymer->chain(1);
    QCOMPARE(chainB->residueCount(), 146);
    QCOMPARE(chainB->sequenceString(), QString("VQLSGEEKAAVLALWDKVNEEEVGGEALGR"
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
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.polymerCount(), 1);
    chemkit::Polymer *polymer = file.polymer();
    QVERIFY(polymer != 0);
    QCOMPARE(polymer->chainCount(), 1);

    chemkit::PolymerChain *chain = polymer->chain(0);
    QCOMPARE(chain->residueCount(), 20);
    QCOMPARE(chain->sequenceString(), QString("ADNRCEQGHILKMFPSTWYV"));
}

void PdbTest::read_fmc()
{
    chemkit::PolymerFile file(dataPath + "fmc.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);
}

QTEST_APPLESS_MAIN(PdbTest)
#include "pdbtest.moc"
