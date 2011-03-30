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

#include <algorithm>

#include <chemkit/polymer.h>
#include <chemkit/polymerfile.h>
#include <chemkit/polymerchain.h>
#include <chemkit/polymerfileformat.h>

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
    std::vector<std::string> formats = chemkit::PolymerFileFormat::formats();
    QVERIFY(std::find(formats.begin(), formats.end(), "pdbml") != formats.end());
}

void PdbmlTest::read_1UBQ()
{
    chemkit::PolymerFile file(dataPath + "1UBQ.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
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
    chemkit::PolymerFile file(dataPath + "2DHB.pdbml");

    bool ok = file.read();
    if(!ok)
        qDebug() << file.errorString();
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
#include "pdbmltest.moc"
