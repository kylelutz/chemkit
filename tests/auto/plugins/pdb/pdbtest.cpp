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

#include <chemkit/chemkit.h>
#include <chemkit/protein.h>
#include <chemkit/nucleicacid.h>
#include <chemkit/proteinchain.h>
#include <chemkit/biochemicalfile.h>
#include <chemkit/nucleicacidchain.h>

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
    QVERIFY(chemkit::BiochemicalFileFormat::formats().contains("pdb"));
}

void PdbTest::read_1BNA()
{
    // create file
    chemkit::BiochemicalFile file;

    // read file
    bool ok = file.read(dataPath + "1BNA.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check nucleic acid
    QCOMPARE(file.nucleicAcidCount(), 1);
    chemkit::NucleicAcid *nucleicAcid = file.nucleicAcid();
    QVERIFY(nucleicAcid != 0);
    QCOMPARE(nucleicAcid->chainCount(), 2);

    // check chains
    chemkit::NucleicAcidChain *chainA = nucleicAcid->chain(0);
    QCOMPARE(chainA->residueCount(), 12);
    QCOMPARE(chainA->sequenceString(), QString("CGCGAATTCGCG"));

    chemkit::NucleicAcidChain *chainB = nucleicAcid->chain(1);
    QCOMPARE(chainB->residueCount(), 12);
    QCOMPARE(chainB->sequenceString(), QString("CGCGAATTCGCG"));
}

void PdbTest::read_1UBQ()
{
    // create file
    chemkit::BiochemicalFile file;

    // read file
    bool ok = file.read(dataPath + "1UBQ.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.proteinCount(), 1);
    chemkit::Protein *protein = file.protein();
    QVERIFY(protein != 0);
    QCOMPARE(protein->chainCount(), 1);

    // check chain
    chemkit::ProteinChain *chain = protein->chain(0);
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
    chemkit::BiochemicalFile file;

    // read file
    bool ok = file.read(dataPath + "2DHB.pdb");
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.proteinCount(), 1);
    chemkit::Protein *protein = file.protein();
    QVERIFY(protein != 0);
    QCOMPARE(protein->chainCount(), 2);

    // chain A
    chemkit::ProteinChain *chainA = protein->chain(0);
    QCOMPARE(chainA->residueCount(), 141);
    QCOMPARE(chainA->sequenceString(), QString("VLSAADKTNVKAAWSKVGGHAGEYGAEALE"
                                               "RMFLGFPTTKTYFPHFDLSHGSAQVKAHGK"
                                               "KVADGLTLAVGHLDDLPGALSDLSNLHAHK"
                                               "LRVDPVNFKLLSHCLLSTLAVHLPNDFTPA"
                                               "VHASLDKFLSSVSTVLTSKYR"));

    // chain B
    chemkit::ProteinChain *chainB = protein->chain(1);
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
    chemkit::BiochemicalFile file(dataPath + "alphabet.pdb");

    // read file
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);

    // check protein
    QCOMPARE(file.proteinCount(), 1);
    chemkit::Protein *protein = file.protein();
    QVERIFY(protein != 0);
    QCOMPARE(protein->chainCount(), 1);
    QCOMPARE(protein->residueCount(), 20);

    chemkit::ProteinChain *chain = protein->chain(0);
    QCOMPARE(chain->sequenceString(), QString("ADNRCEQGHILKMFPSTWYV"));
}

void PdbTest::read_fmc()
{
    chemkit::BiochemicalFile file(dataPath + "fmc.pdb");
    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString();
    QVERIFY(ok);
}

QTEST_APPLESS_MAIN(PdbTest)
#include "pdbtest.moc"
