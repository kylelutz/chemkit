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
#include <QtXml>

#include <algorithm>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/chemicalfile.h>
#include <chemkit/partialchargepredictor.h>

const QString dataPath = "../../../data/";

class MmffTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void validate();
};

void MmffTest::initTestCase()
{
    std::vector<std::string> typers = chemkit::AtomTyper::typers();
    QVERIFY(std::find(typers.begin(), typers.end(), "mmff") != typers.end());

    std::vector<std::string> predictors = chemkit::PartialChargePredictor::predictors();
    QVERIFY(std::find(predictors.begin(), predictors.end(), "mmff") != predictors.end());

    std::vector<std::string> forceFields = chemkit::ForceField::forceFields();
    QVERIFY(std::find(forceFields.begin(), forceFields.end(), "mmff") != forceFields.end());
}

// The validate() method validates the MMFF force field using the MMFF94
// Validation Suite from <http://www.ccl.net/cca/data/MMFF94/>. The suite
// includes 753 molecules and each is check for correct atom typing, atom
// charge assignment, and total energy.
void MmffTest::validate()
{
    // open molecule data file
    chemkit::ChemicalFile dataFile(dataPath + "MMFF94_hypervalent.mol2");
    bool ok = dataFile.read();
    if(!ok)
        qDebug() << dataFile.errorString();
    QVERIFY(ok);
    QCOMPARE(dataFile.moleculeCount(), 753);

    // open expected results file
    QFile expectedFile("mmff94.expected");
    if(!expectedFile.open(QFile::ReadOnly)){
        qDebug() << expectedFile.errorString();
        QFAIL("Failed to open expected data file.");
    }

    QDomDocument expectedFileDocument;
    expectedFileDocument.setContent(&expectedFile);
    QDomElement expectedMolecule = expectedFileDocument.documentElement().firstChildElement();
    QCOMPARE(expectedMolecule.tagName(), QString("molecule"));

    // validate molecules
    QList<chemkit::ForceField *> failedMolecules;
    foreach(const chemkit::Molecule *molecule, dataFile.molecules()){
        bool failed = false;

        // check for correct expected molecule
        QCOMPARE(expectedMolecule.attribute("name").toStdString(), molecule->name());

        // create mmff force field
        chemkit::ForceField *forceField = chemkit::ForceField::create("mmff");
        QVERIFY(forceField);

        // add molecule and setup force field
        forceField->addMolecule(molecule);
        bool setup = forceField->setup();
        if(!setup){
            //failed = true;
        }

        // verify atoms
        int atomCount = forceField->atomCount();
        int expectedAtomCount = expectedMolecule.attribute("atomCount").toInt();
        if(atomCount != expectedAtomCount){
            failed = true;
        }

        QDomElement expectedAtom = expectedMolecule.firstChildElement();
        QCOMPARE(expectedAtom.tagName(), QString("atom"));
        foreach(const chemkit::ForceFieldAtom *forceFieldAtom, forceField->atoms()){
            std::string type = forceFieldAtom->type();
            std::string expectedType = expectedAtom.attribute("type").toStdString();
            if(type != expectedType){
                failed = true;
            }

            double charge = forceFieldAtom->charge();
            double expectedCharge = expectedAtom.attribute("charge").toDouble();
            double chargeDifference = qAbs(charge - expectedCharge);
            if(chargeDifference > 0.1){
                failed = true;
            }

            expectedAtom = expectedAtom.nextSiblingElement();
        }

        // verify energy
        double energy = forceField->energy();
        double expectedEnergy = expectedMolecule.attribute("energy").toDouble();
        double energyDifference = qAbs(energy - expectedEnergy);
        if(energyDifference > 1.0){
            failed = true;
        }

        // move expected molecule to next molecule element
        expectedMolecule = expectedMolecule.nextSiblingElement();

        if(failed){
            failedMolecules.append(forceField);
        }
        else{
            delete forceField;
        }
    }

    // write actual results file if any molecules failed
    if(failedMolecules.size() > 0){
        QFile actualFile("mmff94.actual");
        actualFile.open(QFile::WriteOnly);

        actualFile.write("<molecules>\n");

        foreach(chemkit::ForceField *forceField, failedMolecules){
            const chemkit::Molecule *molecule = forceField->molecules()[0];

            actualFile.write(QString("  <molecule name=\"%1\" energy=\"%2\" atomCount=\"%3\">\n")
                                .arg(molecule->name().c_str())
                                .arg(forceField->energy())
                                .arg(forceField->atomCount())
                                .toAscii());

            foreach(const chemkit::ForceFieldAtom *forceFieldAtom, forceField->atoms()){
                actualFile.write(QString("    <atom type=\"%1\" charge=\"%2\"/>\n")
                                    .arg(forceFieldAtom->type().c_str())
                                    .arg(forceFieldAtom->charge())
                                    .toAscii());
            }

            actualFile.write("  </molecule>\n");
            delete forceField;
        }

        actualFile.write("</molecules>\n");
        actualFile.close();
    }

    // verify that there are no failed molecules
    QCOMPARE(failedMolecules.size(), 0);
}

QTEST_APPLESS_MAIN(MmffTest)
#include "mmfftest.moc"
