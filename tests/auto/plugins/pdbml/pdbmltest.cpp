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
#include <chemkit/biochemicalfile.h>

const QString dataPath = "../../../data/";

class PdbmlTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void read_1UBQ();
        void read_2DHB();
};

void PdbmlTest::initTestCase()
{
    QVERIFY(chemkit::BiochemicalFileFormat::formats().contains("pdbml"));
}

void PdbmlTest::read_1UBQ()
{
    chemkit::BiochemicalFile file(dataPath + "1UBQ.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);
    QCOMPARE(file.proteinCount(), 1);

    // protein
    chemkit::Protein *protein = file.protein();
    QCOMPARE(protein->chainCount(), 1);
    QCOMPARE(protein->residueCount(), 76);

    // chain
    chemkit::ProteinChain *chain = protein->chain(0);
    QVERIFY(chain != 0);
    QCOMPARE(chain->sequenceString(), QString("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"
                                              "DKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
                                              "KESTLHLVLRLRGG"));

    // molecule
    chemkit::Molecule *molecule = protein->molecule();
    QCOMPARE(molecule->atomCount(), 660);
}

void PdbmlTest::read_2DHB()
{
    chemkit::BiochemicalFile file(dataPath + "2DHB.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
    QVERIFY(ok);
    QCOMPARE(file.proteinCount(), 1);

    // protein
    chemkit::Protein *protein = file.protein();
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

QTEST_APPLESS_MAIN(PdbmlTest)
#include "pdbmltest.moc"
