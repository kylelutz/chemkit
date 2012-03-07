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

#include "mmfftest.h"

#include <QtXml>

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>
#include <chemkit/aromaticitymodel.h>
#include <chemkit/moleculardescriptor.h>
#include <chemkit/partialchargepredictor.h>

const std::string dataPath = "../../../data/";

void MmffTest::initTestCase()
{
    // verify that the mmff plugin registered itself correctly
    QVERIFY(boost::count(chemkit::AtomTyper::typers(), "mmff") == 1);
    QVERIFY(boost::count(chemkit::ForceField::forceFields(), "mmff") == 1);
    QVERIFY(boost::count(chemkit::PartialChargePredictor::predictors(), "mmff") == 1);
    QVERIFY(boost::count(chemkit::AromaticityModel::models(), "mmff") == 1);
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "mmff-energy") == 1);
}

// The validate() method validates the MMFF force field using the MMFF94
// Validation Suite from <http://www.ccl.net/cca/data/MMFF94/>. The suite
// includes 753 molecules and each is check for correct atom typing, atom
// charge assignment, and total energy.
void MmffTest::validate()
{
    // open molecule data file
    chemkit::MoleculeFile dataFile(dataPath + "MMFF94_hypervalent.mol2");
    bool ok = dataFile.read();
    if(!ok)
        qDebug() << dataFile.errorString().c_str();
    QVERIFY(ok);
    QCOMPARE(dataFile.moleculeCount(), size_t(753));

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
    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, dataFile.molecules()){
        bool failed = false;

        // check for correct expected molecule
        QByteArray name = expectedMolecule.attribute("name").toAscii();
        QCOMPARE(name.constData(), molecule->name().c_str());

        // create mmff force field
        chemkit::ForceField *forceField = chemkit::ForceField::create("mmff");
        QVERIFY(forceField);

        // add molecule and setup force field
        forceField->setMolecule(molecule.get());
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
            QByteArray expectedType = expectedAtom.attribute("type").toAscii();
            if(type != expectedType.constData()){
                failed = true;
            }

            double charge = forceFieldAtom->charge();
            double expectedCharge = expectedAtom.attribute("charge").toDouble();
            double chargeDifference = std::abs(charge - expectedCharge);
            if(chargeDifference > 0.1){
                failed = true;
            }

            expectedAtom = expectedAtom.nextSiblingElement();
        }

        // verify energy
        double energy = forceField->energy();
        double expectedEnergy = expectedMolecule.attribute("energy").toDouble();
        double energyDifference = std::abs(energy - expectedEnergy);
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
            const chemkit::Molecule *molecule = forceField->molecule();

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
