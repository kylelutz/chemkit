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

#include <chemkit/molecule.h>
#include <chemkit/biochemicalfile.h>
#include <chemkit/proteindatabank.h>

class ProteinDataBankTest : public QObject
{
	Q_OBJECT

	private slots:
        void downloadFile();
        void downloadLigand();
};

void ProteinDataBankTest::downloadFile()
{
    chemkit::ProteinDataBank pdb;

    // download file for lysozyme
    chemkit::BiochemicalFile *file = pdb.downloadFile("2LYZ");
    QVERIFY(file != 0);

    QCOMPARE(file->proteinCount(), 1);
    chemkit::Protein *protein = file->protein();
    QCOMPARE(protein->residueCount(), 129);

    delete file;
}

void ProteinDataBankTest::downloadLigand()
{
    chemkit::ProteinDataBank pdb;

    chemkit::Molecule *molecule = pdb.downloadLigand("ADP");
    QVERIFY(molecule != 0);
    QCOMPARE(molecule->name(), QString("ADP"));
    QCOMPARE(molecule->formula(), QString("C10H15N5O10P2"));
}

QTEST_MAIN(ProteinDataBankTest)
#include "proteindatabanktest.moc"
