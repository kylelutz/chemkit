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

#include "smilestest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/atom.h>
#include <chemkit/ring.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/moleculefile.h>
#include <chemkit/aromaticitymodel.h>
#include <chemkit/substructurequery.h>
#include <chemkit/moleculefileformat.h>

const std::string dataPath = "../../../data/";

void SmilesTest::initTestCase()
{
    // verify that the smiles plugin registered itself correctly
    QVERIFY(boost::count(chemkit::LineFormat::formats(), "smiles") == 1);
    QVERIFY(boost::count(chemkit::AromaticityModel::models(), "daylight") == 1);
    QVERIFY(boost::count(chemkit::MoleculeFileFormat::formats(), "smi") == 1);
}

void SmilesTest::COMPARE_SMILES(const chemkit::Molecule *molecule, const std::string &smiles)
{
    chemkit::Molecule moleculeFromSmiles(smiles, "smiles");

    chemkit::SubstructureQuery query(molecule);
    query.setFlags(chemkit::SubstructureQuery::CompareAromaticity |
                   chemkit::SubstructureQuery::CompareExact);
    bool equal = query.matches(&moleculeFromSmiles);
    if(!equal){
        qDebug() << "Actual SMILES: " << molecule->formula("smiles").c_str();
        qDebug() << "Actual formula: " << moleculeFromSmiles.formula().c_str();
        qDebug() << "Expected SMILES: " << smiles.c_str();
        qDebug() << "Expected formula: " << molecule->formula().c_str();
    }
    QVERIFY(equal);
}

// --- Molecule Tests ------------------------------------------------------ //
void SmilesTest::acenaphthylene()
{
    chemkit::Molecule molecule("c3cc1cccc2C=Cc(c12)c3", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H8"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aceticAcid()
{
    chemkit::Molecule molecule("CC(=O)O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H4O2"));
    QCOMPARE(molecule.bondCount(), size_t(7));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::adenine()
{
    chemkit::Molecule molecule("n1c(c2c(nc1)ncn2)N", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5N5"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::alanine()
{
    chemkit::Molecule molecule("O=C(O)[C@H](N)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C3H7NO2"));
    QCOMPARE(molecule.bondCount(), size_t(12));

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Carbon) && atom->isBondedTo(chemkit::Atom::Nitrogen)){
            QVERIFY(atom->chirality() == chemkit::Stereochemistry::R);
        }
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ampicillin()
{
    chemkit::Molecule molecule("O=C(O)[C@@H]2N3C(=O)[C@@H](NC(=O)[C@@H]"
                               "(c1ccccc1)N)[C@H]3SC2(C)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C16H19N3O4S"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::anthracene()
{
    chemkit::Molecule molecule("c1ccc2cc3ccccc3cc2c1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C14H10"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::anthraquinone()
{
    chemkit::Molecule molecule("O=C2c1ccccc1C(=O)c3ccccc23", "smiles");
    QCOMPARE(molecule.formula(), std::string("C14H8O2"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::arsabenzene()
{
    chemkit::Molecule molecule("[as]1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5As"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::arsole()
{
    chemkit::Molecule molecule("c1[as]ccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H5As"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aspirin()
{
    chemkit::Molecule molecule("O=C(Oc1ccccc1C(=O)O)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C9H8O4"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::aziridine()
{
    chemkit::Molecule molecule("N1CC1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H5N"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1NC1");
}

void SmilesTest::azulene()
{
    chemkit::Molecule molecule("c1cccc2cccc2c1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H8"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzene()
{
    chemkit::Molecule molecule("c1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H6"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1=CC=CC=C1");
}

void SmilesTest::benzofuran()
{
    chemkit::Molecule molecule("o2c1ccccc1cc2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H6O"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzofurazan()
{
    chemkit::Molecule molecule("n1onc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H4N2O"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::benzyne()
{
    chemkit::Molecule molecule("C\\1#C\\C=C/C=C/1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H4"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::binol()
{
    chemkit::Molecule molecule("Oc1c(c2c(O)ccc3c2cccc3)c(cccc4)c4cc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C20H14O2"));
    QCOMPARE(molecule.ringCount(), size_t(4));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), size_t(6));
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biphenyl()
{
    chemkit::Molecule molecule("c1ccccc1(c2ccccc2)", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H10"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biphenylene()
{
    chemkit::Molecule molecule("c3cc2c1c(cccc1)c2cc3", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H8"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->size() == 6){
            QCOMPARE(ring->isAromatic(), true);
        }
        else if(ring->size() == 4){
            QCOMPARE(ring->isAromatic(), false);
        }
    }

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::biperiden()
{
    chemkit::Molecule molecule("OC(c1ccccc1)(CCN2CCCCC2)C4C3\\C=C/C(C3)C4", "smiles");
    QCOMPARE(molecule.formula(), std::string("C21H29NO"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::borinine()
{
    chemkit::Molecule molecule("b1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5B"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::borole()
{
    chemkit::Molecule molecule("C1=CC=CB1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H5B"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::buckminsterfullerene()
{
    chemkit::Molecule molecule("c12c3c4c5c1c6c7c8c2c9c%10c3c%11c%12c4c%13c%14"
                               "c5c%15c6c%16c7c%17c%18c8c9c%19c%20c%10c%11c%21"
                               "c%22c%12c%13c%23c%24c%14c%15c%25c%16c%26c%17"
                               "c%27c%18c%19c%28c%20c%21c%29c%22c%23c%30c%24"
                               "c%25c%26c%31c%27c%28c%29c%30%31", "smiles");
    QCOMPARE(molecule.formula(), std::string("C60"));
}

void SmilesTest::butene()
{
    // cis butene
    chemkit::Molecule cis("C(=C\\C)\\C", "smiles");
    QCOMPARE(cis.formula(), std::string("C4H8"));

    foreach(const chemkit::Bond *bond, cis.bonds()){
        if(bond->order() == chemkit::Bond::Double){
            QVERIFY(bond->stereochemistry() == chemkit::Stereochemistry::E);
        }
    }

    // trans butene
    chemkit::Molecule trans("C(=C/C)\\C", "smiles");
    QCOMPARE(trans.formula(), std::string("C4H8"));

    foreach(const chemkit::Bond *bond, trans.bonds()){
        if(bond->order() == chemkit::Bond::Double){
            QVERIFY(bond->stereochemistry() == chemkit::Stereochemistry::Z);
        }
    }
}

void SmilesTest::caffeine()
{
    chemkit::Molecule molecule("O=C2N(c1ncn(c1C(=O)N2C)C)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H10N4O2"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::camphor()
{
    chemkit::Molecule molecule("O=C1CC2CCC1(C)C2(C)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H16O"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::carbazole()
{
    chemkit::Molecule molecule("c1cccc3c1c2c(cccc2)n3", "smiles");
    QCOMPARE(molecule.ringCount(), size_t(3));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->contains(chemkit::Atom::Nitrogen)){
            QCOMPARE(ring->size(), size_t(5));
        }
        else{
            QCOMPARE(ring->size(), size_t(6));
        }

        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cholesterol()
{
    chemkit::Molecule molecule("O[C@@H]4C/C3=C/C[C@@H]1[C@H](CC[C@]2([C@H]1CC"
                               "[C@@H]2[C@H](C)CCCC(C)C)C)[C@@]3(C)CC4", "smiles");
    QCOMPARE(molecule.formula(), std::string("C27H46O"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::chrysene()
{
    chemkit::Molecule molecule("c4c1c(ccc2ccccc12)c3ccccc3c4", "smiles");
    QCOMPARE(molecule.formula(), std::string("C18H12"));
    QCOMPARE(molecule.ringCount(), size_t(4));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), size_t(6));
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cinnoline()
{
    chemkit::Molecule molecule("n1nccc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H6N2"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::colchicine()
{
    chemkit::Molecule molecule("O=C(N[C@@H]3C\\1=C\\C(=O)C(\\OC)=C/C=C/1c2c"
                               "(cc(OC)c(OC)c2OC)CC3)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C22H25NO6"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::copperSulfate()
{
    chemkit::Molecule molecule("[Cu+2].[O-]S(=O)(=O)[O-]", "smiles");
    QCOMPARE(molecule.formula(), std::string("CuO4S"));
    QCOMPARE(molecule.bondCount(), size_t(4));
    QCOMPARE(molecule.fragmentCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::corannulene()
{
    chemkit::Molecule molecule("c16ccc2ccc3ccc5c4c(c1c2c34)c(cc5)cc6", "smiles");
    QCOMPARE(molecule.formula(), std::string("C20H10"));
    QCOMPARE(molecule.ringCount(), size_t(6));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::coronene()
{
    chemkit::Molecule molecule("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67", "smiles");
    QCOMPARE(molecule.formula(), std::string("C24H12"));
    QCOMPARE(molecule.ringCount(), size_t(7));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cubane()
{
    chemkit::Molecule molecule("C12C3C4C1C5C2C3C45", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H8"));
    QCOMPARE(molecule.bondCount(), size_t(20));
    QCOMPARE(molecule.ringCount(), size_t(5));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cyanide()
{
    chemkit::Molecule molecule("C#N", "smiles");
    QCOMPARE(molecule.formula(), std::string("CHN"));
    QCOMPARE(molecule.bondCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::cytosine()
{
    chemkit::Molecule molecule("O=C1/N=C\\C=C(\\N)N1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H5N3O"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::decalin()
{
    chemkit::Molecule molecule("C1CCC2CCCCC2C1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H18"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dibenzofuran()
{
    chemkit::Molecule molecule("o2c1ccccc1c3c2cccc3", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H8O"));
    QCOMPARE(molecule.ringCount(), size_t(3));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);
    QCOMPARE(molecule.rings()[2]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dichloroethene()
{
    chemkit::Molecule molecule("Cl[C@H]=CCl", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H2Cl2"));
    QCOMPARE(molecule.bondCount(), size_t(5));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dihydrogen()
{
    chemkit::Molecule molecule("[H][H]", "smiles");
    QCOMPARE(molecule.formula(), std::string("H2"));
    QCOMPARE(molecule.bondCount(), size_t(1));

    QCOMPARE(molecule.formula("smiles"), std::string("[H][H]"));
    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::dinitrogen()
{
    chemkit::Molecule molecule("N#N", "smiles");
    QCOMPARE(molecule.formula(), std::string("N2"));

    QCOMPARE(molecule.formula("smiles"), std::string("N#N"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ethane()
{
    chemkit::Molecule molecule("CC", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H6"));

    QCOMPARE(molecule.formula("smiles"), std::string("CC"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::fluorenone()
{
    chemkit::Molecule molecule("O=C3c1ccccc1c2c3cccc2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C13H8O"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::folate()
{
    chemkit::Molecule molecule("O=C(O)[C@@H](NC(=O)c1ccc(cc1)NCc2nc3c"
                               "(nc2)N/C(=N\\C3=O)N)CCC(=O)O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C19H19N7O6"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::furan()
{
    chemkit::Molecule molecule("o1cccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4O"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::furazan()
{
    chemkit::Molecule molecule("n1oncc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H2N2O"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::glucose()
{
    chemkit::Molecule molecule("OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H12O6"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::guanine()
{
    chemkit::Molecule molecule("NC1=Nc2[nH]cnc2C(=O)N1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5N5O"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::heavyWater()
{
    chemkit::Molecule molecule("[2H]O[2H]", "smiles");
    QCOMPARE(molecule.formula(), std::string("H2O"));
    QCOMPARE(molecule.bondCount(), size_t(2));

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(atom->massNumber(), chemkit::Atom::MassNumberType(2));
        }
    }

    QCOMPARE(molecule.formula("smiles"), std::string("[2H]O[2H]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::histidine()
{
    chemkit::Molecule molecule("N[C@@H](Cc1[nH]cnc1)C(O)=O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H9N3O2"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::hydride()
{
    chemkit::Molecule molecule("[H-]", "smiles");
    QCOMPARE(molecule.formula(), std::string("H"));
    QCOMPARE(molecule.bondCount(), size_t(0));
    //QCOMPARE(molecule.atom(0)->formalCharge(), -1);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::hydronium()
{
    chemkit::Molecule molecule("[OH3+]", "smiles");
    QCOMPARE(molecule.formula(), std::string("H3O"));
    QCOMPARE(molecule.bondCount(), size_t(3));

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(atom->formalCharge(), 1);
        }
    }

    QCOMPARE(molecule.formula("smiles"), std::string("[OH3+]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ibuprofen()
{
    chemkit::Molecule molecule("CC(C(=O)O)c1ccc(CC(C)C)cc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C13H18O2"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indazole()
{
    chemkit::Molecule molecule("n2cc1ccccc1n2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C7H6N2"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indene()
{
    chemkit::Molecule molecule("c1cccc2c1\\C=C/C2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C9H8"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indole()
{
    chemkit::Molecule molecule("c1cccc2c1ccn2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H7N"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::indolizine()
{
    chemkit::Molecule molecule("c1ccc2ccccn12", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H7N"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ipratropium()
{
    chemkit::Molecule molecule("O=C(OC2CC1[N+](C)(C(CC1)C2)C(C)C)C(c3ccccc3)CO", "smiles");
    QCOMPARE(molecule.formula(), std::string("C20H30NO3"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isobutane()
{
    chemkit::Molecule molecule("CC(C)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H10"));
    QCOMPARE(molecule.bondCount(), size_t(13));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isoindene()
{
    chemkit::Molecule molecule("C12=CCC=C1C=CC=C2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C9H8"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->size() == 5){
            QCOMPARE(ring->isAromatic(), false);
        }
        else if(ring->size() == 6){
            QCOMPARE(ring->isAromatic(), true);
        }
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::isoindole()
{
    chemkit::Molecule molecule("c1cccc2c1cnc2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H7N"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::melatonin()
{
    chemkit::Molecule molecule("O=C(NCCc2c1cc(OC)ccc1nc2)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C13H16N2O2"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::naphthalene()
{
    chemkit::Molecule molecule("c1ccc2ccccc2c1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H8"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::nicotine()
{
    chemkit::Molecule molecule("CN1CCC[C@H]1c2cccnc2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H14N2"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::nitrobenzene()
{
    chemkit::Molecule molecule("[O-][N+](=O)c1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H5NO2"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::ovalene()
{
    chemkit::Molecule molecule("c1cc2c3c4c1ccc5cc6c7c8c(ccc9=c8c1c(cc9)cc"
                               "(c3c1c7c54)cc2)cc6", "smiles");
    QCOMPARE(molecule.formula(), std::string("C32H14"));
    QCOMPARE(molecule.ringCount(), size_t(10));
}

void SmilesTest::oxazole()
{
    chemkit::Molecule molecule("n1ccoc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C3H3NO"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pentacene()
{
    chemkit::Molecule molecule("c45cc3cc2cc1ccccc1cc2cc3cc4cccc5", "smiles");
    QCOMPARE(molecule.formula(), std::string("C22H14"));
    QCOMPARE(molecule.ringCount(), size_t(5));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), size_t(6));
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pentalene()
{
    chemkit::Molecule molecule("c1cc2cccc2c1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H6"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::perylene()
{
    chemkit::Molecule molecule("c1ccc5cccc4c5c1c2cccc3cccc4c23", "smiles");
    QCOMPARE(molecule.formula(), std::string("C20H12"));
    QCOMPARE(molecule.ringCount(), size_t(5));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenanthrene()
{
    chemkit::Molecule molecule("c1ccc2c(c1)ccc3ccccc32", "smiles");
    QCOMPARE(molecule.formula(), std::string("C14H10"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), size_t(6));
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenothiazine()
{
    chemkit::Molecule molecule("c1ccc2Sc3ccccc3Nc2c1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H9NS"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phenoxazine()
{
    chemkit::Molecule molecule("O2c1ccccc1Nc3c2cccc3", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H9NO"));
    QCOMPARE(molecule.ringCount(), size_t(3));

    //COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phosphole()
{
    chemkit::Molecule molecule("c1cccp1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H5P"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phosphorine()
{
    chemkit::Molecule molecule("p1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5P"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::phthalimide()
{
    chemkit::Molecule molecule("O=C2c1ccccc1C(=O)N2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H5NO2"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::porphin()
{
    chemkit::Molecule molecule("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3", "smiles");
    QCOMPARE(molecule.formula(), std::string("C20H14N4"));
    QCOMPARE(molecule.ringCount(), size_t(5));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::proline()
{
    chemkit::Molecule molecule("O=C(O)C1NCCC1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H9NO2"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::proton()
{
    chemkit::Molecule molecule("[H+]", "smiles");
    QCOMPARE(molecule.formula(), std::string("H"));
    QCOMPARE(molecule.bondCount(), size_t(0));

    QCOMPARE(molecule.formula("smiles"), std::string("[H+]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::purine()
{
    chemkit::Molecule molecule("n1cc2c(nc1)ncn2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H4N4"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyranium()
{
    chemkit::Molecule molecule("[o+]1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5O"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    foreach(const chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(atom->formalCharge(), 1);
        }
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrazole()
{
    chemkit::Molecule molecule("n1cccn1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C3H4N2"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrene()
{
    chemkit::Molecule molecule("c3ccc2ccc1cccc4c1c2c3cc4", "smiles");
    QCOMPARE(molecule.formula(), std::string("C16H10"));
    QCOMPARE(molecule.ringCount(), size_t(4));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyridazine()
{
    chemkit::Molecule molecule("n1ncccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4N2"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyridine()
{
    chemkit::Molecule molecule("n1ccccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H5N"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrimidine()
{
    chemkit::Molecule molecule("n1cccnc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4N2"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::pyrrole()
{
    chemkit::Molecule molecule("n1cccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H5N"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::quinoline()
{
    chemkit::Molecule molecule("n1cccc2ccccc12", "smiles");
    QCOMPARE(molecule.formula(), std::string("C9H7N"));
    QCOMPARE(molecule.ringCount(), size_t(2));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);
    QCOMPARE(molecule.rings()[1]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::rhodizonicAcid()
{
    chemkit::Molecule molecule("O=C1C(/O)=C(/O)C(=O)C(=O)C1=O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C6H2O6"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::selenophene()
{
    chemkit::Molecule molecule("[se]1cccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4Se"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::sodiumChloride()
{
    chemkit::Molecule molecule("[Na+].[Cl-]", "smiles");
    QCOMPARE(molecule.formula(), std::string("ClNa"));
    QCOMPARE(molecule.bondCount(), size_t(0));
    QCOMPARE(molecule.fragmentCount(), size_t(2));

    QCOMPARE(molecule.formula("smiles"), std::string("[Na+].[Cl-]"));
    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::stilbene()
{
    chemkit::Molecule molecule("c2(\\C=C\\c1ccccc1)ccccc2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C14H12"));
    QCOMPARE(molecule.bondCount(), size_t(27));
    QCOMPARE(molecule.ringCount(), size_t(2));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        QCOMPARE(ring->size(), size_t(6));
        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::sulfurHexafluoride()
{
    chemkit::Molecule molecule("FS(F)(F)(F)(F)F", "smiles");
    QCOMPARE(molecule.formula(), std::string("F6S"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::taxol()
{
    chemkit::Molecule molecule("O=C(N[C@@H](c1ccccc1)[C@@H](O)C(=O)O[C@H]5C"
                               "[C@@]6(O)[C@@H](OC(=O)c2ccccc2)[C@H]3[C@@](C)"
                               "([C@@H](O)C[C@H]4OC[C@@]34OC(C)=O)C(=O)[C@H]"
                               "(OC(C)=O)\\C(=C5/C)[C@]6(C)C)c7ccccc7", "smiles");
    QCOMPARE(molecule.formula(), std::string("C47H51NO14"));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::tetraphenylene()
{
    chemkit::Molecule molecule("c5cc4c1c(cccc1)c2ccccc2c3ccccc3c4cc5", "smiles");
    QCOMPARE(molecule.formula(), std::string("C24H16"));
    QCOMPARE(molecule.ringCount(), size_t(5));
}

void SmilesTest::tetralin()
{
    chemkit::Molecule molecule("c1ccc2c(c1)CCCC2", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H12"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiamin()
{
    chemkit::Molecule molecule("n1c(c(cnc1C)C[n+]2c(c(sc2)CCO)C)N", "smiles");
    QCOMPARE(molecule.formula(), std::string("C12H16N4OS"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiirane()
{
    chemkit::Molecule molecule("C1CS1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H4S"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thiophene()
{
    chemkit::Molecule molecule("s1cccc1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4S"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thujone()
{
    chemkit::Molecule molecule("C[C@@H]([C@@H](C2)[C@]2([C@@H](C)C)C1)C1=O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C10H16O"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::thymine()
{
    chemkit::Molecule molecule("O=C1\\C(=C/NC(=O)N1)C", "smiles");
    QCOMPARE(molecule.formula(), std::string("C5H6N2O2"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::triazole()
{
    chemkit::Molecule molecule("n1ccnn1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C2H3N3"));
    QCOMPARE(molecule.ringCount(), size_t(1));
    QCOMPARE(molecule.rings()[0]->isAromatic(), true);

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::triphenylene()
{
    chemkit::Molecule molecule("c4cc3c1c(cccc1)c2ccccc2c3cc4", "smiles");
    QCOMPARE(molecule.formula(), std::string("C18H12"));
    QCOMPARE(molecule.ringCount(), size_t(4));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
    COMPARE_SMILES(&molecule, "C1(C=CC=C3)=C3C(C=CC=C4)=C4C2=C1C=CC=C2");
}

void SmilesTest::tropone()
{
    chemkit::Molecule molecule("C1=CC=CC(=O)C=C1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C7H6O"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::tryptophan()
{
    chemkit::Molecule molecule("N[C@@H](Cc1c2ccccc2nc1)C(O)=O", "smiles");
    QCOMPARE(molecule.formula(), std::string("C11H12N2O2"));
    QCOMPARE(molecule.ringCount(), size_t(2));

    foreach(const chemkit::Ring *ring, molecule.rings()){
        if(ring->contains(chemkit::Atom::Nitrogen)){
            QCOMPARE(ring->size(), size_t(5));
        }
        else{
            QCOMPARE(ring->size(), size_t(6));
        }

        QCOMPARE(ring->isAromatic(), true);
    }

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::uracil()
{
    chemkit::Molecule molecule("O=C1\\C=C/NC(=O)N1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C4H4N2O2"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

void SmilesTest::vanillin()
{
    chemkit::Molecule molecule("O=CC1=CC(OC)=C(O)C=C1", "smiles");
    QCOMPARE(molecule.formula(), std::string("C8H8O3"));
    QCOMPARE(molecule.ringCount(), size_t(1));

    COMPARE_SMILES(&molecule, molecule.formula("smiles"));
}

// --- Feature Tests ------------------------------------------------------- //
void SmilesTest::addHydrogens()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    // default is true
    QCOMPARE(format->option("add-implicit-hydrogens").toBool(), true);
    chemkit::Molecule *molecule = format->read("C");
    QCOMPARE(molecule->formula(), std::string("CH4"));
    delete molecule;

    format->setOption("add-implicit-hydrogens", false);
    molecule = format->read("C");
    QCOMPARE(molecule->formula(), std::string("C"));
    delete molecule;

    delete format;
}

void SmilesTest::isotope()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("[14CH4]");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), std::string("CH4"));

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(atom->massNumber(), chemkit::Atom::MassNumberType(14));
        }
    }

    delete molecule;

    molecule = format->read("[238U]");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), std::string("U"));
    QCOMPARE(molecule->atom(0)->massNumber(), chemkit::Atom::MassNumberType(238));

    delete molecule;
    delete format;
}

void SmilesTest::kekulize()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    // default is false
    QCOMPARE(format->option("kekulize").toBool(), false);

    chemkit::Molecule benzene("c1ccccc1", "smiles");
    QCOMPARE(format->write(&benzene), std::string("c1ccccc1"));

    format->setOption("kekulize", true);
    QCOMPARE(format->write(&benzene), std::string("C1=CC=CC=C1"));

    delete format;
}

void SmilesTest::quadrupleBond()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C$C");
    QVERIFY(molecule);
    QCOMPARE(molecule->formula(), std::string("C2"));
    QCOMPARE(molecule->bondCount(), size_t(1));
    QCOMPARE(molecule->bonds()[0]->order(), chemkit::Bond::BondOrderType(4));

    delete molecule;
    delete format;
}

// --- Invalid Tests ------------------------------------------------------- //
void SmilesTest::extraParenthesis()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C(C=O))C");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().empty() == false);

    delete format;
}

void SmilesTest::invalidAtom()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("CCX");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().empty() == false);

    delete format;
}

void SmilesTest::wildcardAtom()
{
    chemkit::LineFormat *format = chemkit::LineFormat::create("smiles");
    QVERIFY(format);

    chemkit::Molecule *molecule = format->read("C*C");
    QVERIFY(molecule == 0);
    QVERIFY(format->errorString().empty() == false);

    delete format;
}

// --- File Tests ---------------------------------------------------------- //
void SmilesTest::herg()
{
    chemkit::MoleculeFile file(dataPath + "herg.smi");

    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(31));
    QCOMPARE(file.molecule(0)->name(), std::string("Amitriptyline"));
    QCOMPARE(file.molecule(0)->formula(), std::string("C20H23N"));
    QCOMPARE(file.molecule(30)->name(), std::string ("Verapamil"));
    QCOMPARE(file.molecule(30)->formula(), std::string("C27H38N2O4"));
}

void SmilesTest::cox2()
{
    chemkit::MoleculeFile file(dataPath + "cox2.smi");

    bool ok = file.read();
    if(!ok)
        qDebug() << "Failed to read file: " << file.errorString().c_str();
    QVERIFY(ok);

    QCOMPARE(file.moleculeCount(), size_t(128));
    QCOMPARE(file.molecule(0)->formula(), std::string("C13H18N2O5S"));
    QCOMPARE(file.molecule(2)->formula(), std::string("C16H13F2NO3S2"));
    QCOMPARE(file.molecule(127)->formula(), std::string("C21H19NO5S"));
}

QTEST_APPLESS_MAIN(SmilesTest)
