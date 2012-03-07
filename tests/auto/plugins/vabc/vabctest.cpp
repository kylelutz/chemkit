/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "vabctest.h"

#include <boost/range/algorithm.hpp>

#include <chemkit/molecule.h>
#include <chemkit/moleculardescriptor.h>

void VabcTest::initTestCase()
{
    // verify that the vabc plugin registered itself correctly
    QVERIFY(boost::count(chemkit::MolecularDescriptor::descriptors(), "vabc") == 1);
}

void VabcTest::vabc_data()
{
    QTest::addColumn<QString>("smilesString");
    QTest::addColumn<QString>("formulaString");
    QTest::addColumn<double>("vabc");

    // general compounds
    QTest::newRow("propionicacid") << "CCC(O)=O" << "C3H6O2" << 75.39;
    QTest::newRow("acetic acid") << "CC(O)=O" << "C2H4O2" << 58.09;
    QTest::newRow("butyricacid") << "CCCC(O)=O" << "C4H8O2" << 92.68;
    QTest::newRow("butanol") << "CCCCO" << "C4H10O" << 86.53;
    QTest::newRow("heptanol") << "CCCCCCCO" << "C7H16O" << 138.42;
    QTest::newRow("octanol") << "CCCCCCCCO" << "C8H18O" << 155.71;
    QTest::newRow("pentanol") << "CCCCCO" << "C5H12O" << 103.83;
    QTest::newRow("s-butanol") << "CCC(C)O" << "C4H10O" << 86.53;
    QTest::newRow("isopropyl alcohol") << "CC(C)O" << "C3H8O" << 69.23;
    QTest::newRow("cyclohexanol") << "OC1CCCCC1" << "C6H12O" << 108.77;
    QTest::newRow("ethanol") << "CCO" << "C2H6O" << 51.94;
    QTest::newRow("ethane-1,2-diol") << "OCCO" << "C2H6O2" << 60.73;
    QTest::newRow("methanol") << "CO" << "CH4O" << 34.64;
    QTest::newRow("allylalcohol") << "OCC=C" << "C3H6O" << 66.60;
    QTest::newRow("propanol") << "CCCO" << "C3H8O" << 69.23;
    QTest::newRow("butane") << "CCCC" << "C4H10" << 77.74;
    QTest::newRow("cyclohexane") << "C1CCCCC1" << "C6H12" << 99.98;
    QTest::newRow("cyclooctane") << "C1CCCCCCC1" << "C8H16" << 134.57;
    QTest::newRow("cyclopentane") << "C1CCCC1" << "C5H10" << 82.68;
    QTest::newRow("cyclopropane") << "C1CC1" << "C3H6" << 48.09;
    QTest::newRow("ethane") << "CC" << "C2H6" << 43.15;
    QTest::newRow("hexane") << "CCCCCC" << "C6H14" << 112.33;
    QTest::newRow("methylcyclopentane") << "CC1CCCC1" << "C6H12" << 99.98;
    QTest::newRow("methane") << "C" << "CH4" << 25.85;
    QTest::newRow("isopentane") << "CCC(C)C" << "C5H12" << 95.04;
    QTest::newRow("octane") << "CCCCCCCC" << "C8H18" << 146.92;
    QTest::newRow("pentane") << "CCCCC" << "C5H12" << 95.04;
    QTest::newRow("propane") << "CCC" << "C3H8" << 60.44;
    QTest::newRow("2-butene-cis") << "C\\C=C/C" << "C4H8" << 75.10;
    QTest::newRow("trimethylethylene") << "CC=C(C)C" << "C5H10" << 92.40;
    QTest::newRow("1,3-butadiene") << "C=CC=C" << "C4H6" << 72.47;
    QTest::newRow("1-butene") << "CCC=C" << "C4H8" << 75.10;
    QTest::newRow("cyclohexene") << "C1CCC=CC1" << "C6H10" << 97.34;
    QTest::newRow("cyclopentene") << "C1CC=CC1" << "C5H8" << 80.04;
    QTest::newRow("chloroethylene") << "ClC=C" << "C2H3Cl" << 55.72;
    QTest::newRow("1-hexene") << "CCCCC=C" << "C6H12" << 109.70;
    QTest::newRow("1-methylcyclohexene") << "CC1=CCCCC1" << "C7H12" << 114.64;
    QTest::newRow("1-octene") << "CCCCCCC=C" << "C8H16" << 144.29;
    QTest::newRow("1,4-pentadiene") << "C=CCC=C" << "C5H8" << 89.76;
    QTest::newRow("1-pentene") << "CCCC=C" << "C5H10" << 92.40;
    QTest::newRow("allene") << "C=C=C" << "C3H4" << 55.17;
    QTest::newRow("propylene") << "CC=C" << "C3H6" << 57.81;
    QTest::newRow("2-butyne") << "CC#CC" << "C4H6" << 72.47;
    QTest::newRow("1-butyne") << "CCC#C" << "C4H6" << 72.47;
    QTest::newRow("acetylene") << "C#C" << "C2H2" << 37.88;
    QTest::newRow("1-hexyne") << "CCCCC#C" << "C6H10" << 107.06;
    QTest::newRow("1-octyne") << "CCCCCCC#C" << "C8H14" << 141.65;
    QTest::newRow("1-pentyne") << "CCCC#C" << "C5H8" << 89.76;
    QTest::newRow("methylacetylene") << "CC#C" << "C3H4" << 55.17;
    QTest::newRow("acetamide") << "CC(=O)N" << "C2H5NO" << 60.30;
    QTest::newRow("n,n-dimethylacetamide") << "CN(C)C(=O)C" << "C4H9NO" << 94.89;
    QTest::newRow("dimethylformamide") << "CN(C)C=O" << "C3H7NO" << 77.59;
    QTest::newRow("formamide") << "NC=O" << "CH3NO" << 43.00;
    QTest::newRow("n-methylacetamide") << "CNC(=O)C" << "C3H7NO" << 77.59;
    QTest::newRow("dibutylamine") << "CCCCNCCCC" << "C8H19N" << 157.92;
    QTest::newRow("1-butylamine") << "CCCCN" << "C4H11N" << 88.74;
    QTest::newRow("ethylene diamine") << "NCCN" << "C2H8N2" << 65.14;
    QTest::newRow("diethylamine") << "CCNCC" << "C4H11N" << 88.74;
    QTest::newRow("triethylamine") << "CCN(CC)CC" << "C6H15N" << 123.33;
    QTest::newRow("ethylamine") << "CCN" << "C2H7N" << 54.15;
    QTest::newRow("hexylamine") << "CCCCCCN" << "C6H15N" << 123.33;
    QTest::newRow("i-propylamine") << "CC(C)N" << "C3H9N" << 71.44;
    QTest::newRow("dimethylamine") << "CNC" << "C2H7N" << 54.15;
    QTest::newRow("trimethylamine") << "CN(C)C" << "C3H9N" << 71.44;
    QTest::newRow("methylamine") << "CN" << "CH5N" << 36.85;
    QTest::newRow("n-methylpiperidine") << "CN1CCCCC1" << "C6H13N" << 110.97;
    QTest::newRow("piperazine") << "C1CNCCN1" << "C4H10N2" << 87.38;
    QTest::newRow("piperidine") << "C1CCNCC1" << "C5H11N" << 93.68;
    QTest::newRow("dipropylamine") << "CCCNCCC" << "C6H15N" << 123.33;
    QTest::newRow("propylamine") << "CCCN" << "C3H9N" << 71.44;
    QTest::newRow("pyrrolidine") << "C1CCNC1" << "C4H9N" << 76.38;
    QTest::newRow("acetophenone") << "CC(=O)c1ccccc1" << "C8H8O" << 121.91;
    QTest::newRow("anthracene") << "c3ccc2cc1ccccc1cc2c3" << "C14H10" << 162.48;
    QTest::newRow("benzene") << "c1ccccc1" << "C6H6" << 81.17;
    QTest::newRow("benzophenone") << "O=C(c1ccccc1)c2ccccc2" << "C13H10O" << 177.23;
    QTest::newRow("benzylalcohol") << "OCc1ccccc1" << "C7H8O" << 107.25;
    QTest::newRow("biphenyl") << "c1ccc(cc1)c2ccccc2" << "C12H10" << 153.78;
    QTest::newRow("ethylbenzene") << "CCc1ccccc1" << "C8H10" << 115.76;
    QTest::newRow("furan") << "c1ccoc1" << "C4H4O" << 58.00;
    QTest::newRow("m-hydroxybenzaldehyde") << "Oc1cccc(C=O)c1" << "C7H6O2" << 113.41;
    QTest::newRow("m-xylene") << "Cc1cccc(C)c1" << "C8H10" << 115.76;
    QTest::newRow("naphthalene") << "c2ccc1ccccc1c2" << "C10H8" << 121.82;
    QTest::newRow("o-methylphenol") << "Cc1ccccc1O" << "C7H8O" << 107.25;
    QTest::newRow("o-xylene") << "Cc1ccccc1C" << "C8H10" << 115.76;
    QTest::newRow("p-methylphenol") << "Cc1ccc(O)cc1" << "C7H8O" << 107.25;
    QTest::newRow("p-hydroxybenzaldehyde") << "Oc1ccc(C=O)cc1" << "C7H6O2" << 113.41;
    QTest::newRow("p-xylene") << "Cc1ccc(C)cc1" << "C8H10" << 115.76;
    QTest::newRow("diphenylmethane") << "C(c1ccccc1)c2ccccc2" << "C13H12" << 171.07;
    QTest::newRow("benzaldehyde") << "O=Cc1ccccc1" << "C7H6O" << 104.62;
    QTest::newRow("benzoicacid") << "OC(=O)c1ccccc1" << "C7H6O2" << 113.41;
    QTest::newRow("phenol") << "Oc1ccccc1" << "C6H6O" << 89.96;
    QTest::newRow("styrene") << "C=Cc1ccccc1" << "C8H8" << 113.12;
    QTest::newRow("p-t-butylphenol") << "CC(C)(C)c1ccc(O)cc1" << "C10H14O" << 159.14;
    QTest::newRow("toluene") << "Cc1ccccc1" << "C7H8" << 98.46;
    QTest::newRow("aniline") << "Nc1ccccc1" << "C6H7N" << 92.16;
    QTest::newRow("benzylamine") << "NCc1ccccc1" << "C7H9N" << 109.46;
    QTest::newRow("quinoline") << "c2ccc1ncccc1c2" << "C9H7N" << 115.52;
    QTest::newRow("imidazole") << "c1c[nH]cn1" << "C3H4N2" << 53.91;
    QTest::newRow("3-methylpyridine") << "Cc1cccnc1" << "C6H7N" << 92.16;
    QTest::newRow("2,6-lutidine") << "Cc1cccc(C)n1" << "C7H9N" << 109.46;
    QTest::newRow("2-methylpyrazine") << "Cc1cnccn1" << "C5H6N2" << 85.86;
    QTest::newRow("3-methylindole") << "Cc1c[nH]c2ccccc12" << "C9H9N" << 118.16;
    QTest::newRow("2-picoline") << "Cc1ccccn1" << "C6H7N" << 92.16;
    QTest::newRow("4-methylpyridine") << "Cc1ccncc1" << "C6H7N" << 92.16;
    QTest::newRow("pyrazine") << "c1cnccn1" << "C4H4N2" << 68.57;
    QTest::newRow("pyridazine") << "c1ccnnc1" << "C4H4N2" << 68.57;
    QTest::newRow("pyridine") << "c1ccncc1" << "C5H5N" << 74.87;
    QTest::newRow("pyrimidine") << "c1cncnc1" << "C4H4N2" << 68.57;
    QTest::newRow("pyrrole") << "c1cc[nH]c1" << "C4H5N" << 60.21;
    QTest::newRow("2-heptanone") << "CCCCCC(C)=O" << "C7H14O" << 135.78;
    QTest::newRow("2-octanone") << "CCCCCCC(C)=O" << "C8H16O" << 153.08;
    QTest::newRow("2-pentanone") << "CCCC(C)=O" << "C5H10O" << 101.19;
    QTest::newRow("3-pentanone") << "CCC(=O)CC" << "C5H10O" << 101.19;
    QTest::newRow("acetaldehyde") << "CC=O" << "C2H4O" << 49.30;
    QTest::newRow("propenal") << "C=CC=O" << "C3H4O" << 63.96;
    QTest::newRow("butyraldehyde") << "CCCC=O" << "C4H8O" << 83.89;
    QTest::newRow("2-butanone") << "CCC(C)=O" << "C4H8O" << 83.89;
    QTest::newRow("formaldehyde") << "C=O" << "CH2O" << 32.01;
//    QTest::newRow("quinone") << "O=C1C=CC(=O)C=C1" << "C6H4O2" << 107.01;
//    QTest::newRow("cyclohexanone") << "O=C1CCCCC1" << "C6H10O" << 114.92;
    QTest::newRow("hexaldehyde") << "CCCCCC=O" << "C6H12O" << 118.49;
    QTest::newRow("2-butanone,3-methyl") << "CC(C)C(C)=O" << "C5H10O" << 101.19;
    QTest::newRow("octanal") << "CCCCCCCC=O" << "C8H16O" << 153.08;
    QTest::newRow("propionaldehyde") << "CCC=O" << "C3H6O" << 66.60;
    QTest::newRow("acetone") << "CC(C)=O" << "C3H6O" << 66.60;
    QTest::newRow("ethylformate") << "CCOC=O" << "C3H6O2" << 75.39;
    QTest::newRow("methylformate") << "COC=O" << "C2H4O2" << 58.09;
    QTest::newRow("dioxane") << "C1COCCO1" << "C4H8O2" << 82.96;
    QTest::newRow("1,3-dioxolane") << "C1COCO1" << "C3H6O2" << 65.67;
    QTest::newRow("ethylether") << "CCOCC" << "C4H10O" << 86.53;
    QTest::newRow("di-i-propylether") << "CC(C)OC(C)C" << "C6H14O" << 121.12;
    QTest::newRow("methylpropylether") << "CCCOC" << "C4H10O" << 86.53;
    QTest::newRow("methyl-t-butylether") << "COC(C)(C)C" << "C5H12O" << 103.83;
    QTest::newRow("dimethylether") << "COC" << "C2H6O" << 51.94;
    QTest::newRow("methoxyethanol") << "COCCO" << "C3H8O2" << 78.02;
    QTest::newRow("dipropylether") << "CCCOCCC" << "C6H14O" << 121.12;
    QTest::newRow("tetrahydrofuran") << "C1CCOC1" << "C4H8O" << 74.17;
    QTest::newRow("tetrahydropyran") << "C1CCOCC1" << "C5H10O" << 91.47;
    QTest::newRow("4-bromo-chlorobenzene") << "c1cc(Br)ccc1Cl" << "C6H4BrCl" << 115.66;
    QTest::newRow("ethylbromide") << "CCBr" << "C2H5Br" << 62.43;
    QTest::newRow("1-bromoheptane") << "CCCCCCCBr" << "C7H15Br" << 148.91;
    QTest::newRow("1-bromohexane") << "CCCCCCBr" << "C6H13Br" << 131.62;
    QTest::newRow("2-bromo-2-methylpropane") << "CC(C)(Br)C" << "C4H9Br" << 97.02;
    QTest::newRow("o-bromonitrobenzene") << "Brc1ccccc1N(=O)=O" << "C6H4BrNO2" << 126.39;
    QTest::newRow("1-bromooctane") << "CCCCCCCCBr" << "C8H17Br" << 166.21;
    QTest::newRow("1-bromopentane") << "CCCCCBr" << "C5H11Br" << 114.32;
    QTest::newRow("2-promopropane") << "CC(C)Br" << "C3H7Br" << 79.73;
    QTest::newRow("p-bromophenol") << "Oc1ccc(Br)cc1" << "C6H5BrO" << 109.24;
//    QTest::newRow("1-chlorobutane") << "CCCCCl" << "C4H9Cl" << 89.00;
    QTest::newRow("2,2,2-trifluoroethanol") << "OCC(F)(F)F" << "C2H3F3O" << 70.14;
    QTest::newRow("dibromomethane") << "BrCBr" << "CH2Br2" << 64.42;
    QTest::newRow("methylenechloride") << "ClCCl" << "CH2Cl2" << 56.27;
    QTest::newRow("chloroiodomethane") << "ICCl" << "CH2ClI" << 66.34;
    QTest::newRow("methylbromide") << "CBr" << "CH3Br" << 45.14;
    QTest::newRow("1-iodobutane") << "CCCCI" << "C4H9I" << 103.02;
    QTest::newRow("1-iodopropane") << "CCCI" << "C3H7I" << 85.72;
    QTest::newRow("methylchloride") << "CCl" << "CH3Cl" << 41.06;
    QTest::newRow("methyliodide") << "CI" << "CH3I" << 51.13;
    QTest::newRow("cis-1,2-dichlorethylene") << "ClC=CCl" << "C2H2Cl2" << 70.93;
    QTest::newRow("bromochloromethane") << "ClCBr" << "CH2BrCl" << 60.35;
    QTest::newRow("chlorodifluoromethane") << "FC(F)Cl" << "CHClF2" << 53.20;
//    QTest::newRow("1-chloro-1,1-difluoroethane") << "CCCl(F)(F)" << "C2H3ClF2" << 70.49;
    QTest::newRow("chlorofluoromethane") << "FCCl" << "CH2ClF" << 47.13;
    QTest::newRow("1,4-dibromobenzene") << "c1cc(Br)ccc1Br" << "C6H4Br2" << 119.73;
    QTest::newRow("1,2-dibromoethane") << "BrCCBr" << "C2H4Br2" << 81.72;
    QTest::newRow("1,2-dibromo-ethylene, 1,2-cis") << "BrC=CBr" << "C2H2Br2" << 79.08;
    QTest::newRow("1,2-dibromopropane") << "CC(Br)CBr" << "C3H6Br2" << 99.01;
    QTest::newRow("dichlorofluoromethane") << "FC(Cl)Cl" << "CHCl2F" << 62.34;
    QTest::newRow("diiodomethane") << "ICI" << "CH2I2" << 76.41;
    QTest::newRow("ethylchloride") << "CCCl" << "C2H5Cl" << 58.36;
    QTest::newRow("1,1-dichloroethane") << "CC(Cl)Cl" << "C2H4Cl2" << 73.57;
    QTest::newRow("fluorobenzene") << "Fc1ccccc1" << "C6H5F" << 87.23;
    QTest::newRow("methylfluoride") << "CF" << "CH3F" << 31.92;
    QTest::newRow("1,1-dichloroethylene") << "ClC(Cl)=C" << "C2H2Cl2" << 70.93;
    QTest::newRow("iodobenzene") << "Ic1ccccc1" << "C6H5I" << 106.44;
    QTest::newRow("ethyliodide") << "CCI" << "C2H5I" << 68.43;
    QTest::newRow("2-iodopropane") << "CC(C)I" << "C3H7I" << 85.72;
    QTest::newRow("2-chloropropane") << "CC(C)Cl" << "C3H7Cl" << 75.66;
    QTest::newRow("1,3-dichlorobenzene") << "Clc1cccc(Cl)c1" << "C6H4Cl2" << 111.59;
    QTest::newRow("3-bromotoluene") << "c1ccc(Br)cc1C" << "C7H7Br" << 117.75;
    QTest::newRow("m-dibromobenzene") << "Brc1cccc(Br)c1" << "C6H4Br2" << 119.73;
    QTest::newRow("o-dichlorobenzene") << "Clc1ccccc1Cl" << "C6H4Cl2" << 111.59;
    QTest::newRow("p-dichlorobenzene") << "Clc1ccc(Cl)cc1" << "C6H4Cl2" << 111.59;
    QTest::newRow("1-chloropentane") << "CCCCCCl" << "C5H11Cl" << 110.25;
    QTest::newRow("chlorobenzene") << "Clc1ccccc1" << "C6H5Cl" << 96.38;
    QTest::newRow("1-chloropropane") << "CCCCl" << "C3H7Cl" << 75.66;
    QTest::newRow("1,2,4,5-tetrafluorobenzene") << "Fc1c(F)cc(F)c(F)c1" << "C6H2F4" << 105.44;
    QTest::newRow("1,2,3,5-tetrafluorobenzene") << "Fc1c(F)c(F)cc(F)c1" << "C6H2F4" << 105.44;
    QTest::newRow("trans-1,2-dichloroethene") << "ClC=CCl" << "C2H2Cl2" << 70.93;
    QTest::newRow("1,3,5-tribromobenzene") << "Brc1cc(Br)cc(Br)c1" << "C6H3Br3" << 139.02;
    QTest::newRow("m-cyanophenol") << "Oc1cccc(C#N)c1" << "C7H5NO" << 112.98;
    QTest::newRow("acetonitrile") << "CC#N" << "C2H3N" << 48.87;
    QTest::newRow("p-cyanophenol") << "Oc1ccc(C#N)cc1" << "C7H5NO" << 112.98;
    QTest::newRow("benzonitrile") << "N#Cc1ccccc1" << "C7H5N" << 104.19;
    QTest::newRow("acrylonitrile") << "C=CC#N" << "C3H3N" << 63.53;
    QTest::newRow("benzene,3-cyano-1-nitro") << "O=N(=O)c1cccc(C#N)c1" << "C7H4N2O2" << 130.13;
    QTest::newRow("nitroethane") << "CCN(=O)=O" << "C2H5NO2" << 69.09;
    QTest::newRow("2-nitropropane") << "CC(C)N(=O)=O" << "C3H7NO2" << 86.39;
    QTest::newRow("m-nitrophenol") << "Oc1cccc(c1)N(=O)=O" << "C6H5NO3" << 115.90;
    QTest::newRow("m-nitrotoluene") << "Cc1cccc(c1)N(=O)=O" << "C7H7NO2" << 124.40;
    QTest::newRow("nitromethane") << "CN(=O)=O" << "CH3NO2" << 51.79;
    QTest::newRow("o-nitrotoluene") << "Cc1ccccc1N(=O)=O" << "C7H7NO2" << 124.40;
    QTest::newRow("p-nitrophenol") << "Oc1ccc(cc1)N(=O)=O" << "C6H5NO3" << 115.90;
    QTest::newRow("nitrobenzene") << "O=N(=O)c1ccccc1" << "C6H5NO2" << 107.11;
    QTest::newRow("1-nitropropane") << "CCCN(=O)=O" << "C3H7NO2" << 86.39;
    QTest::newRow("carbonmonoxide") << "[C-]#[O+]" << "CO" << 29.37;
    QTest::newRow("carbondioxide") << "O=C=O" << "CO2" << 38.16;
    QTest::newRow("hydrogen") << "[H][H]" << "H2" << 8.56;
    QTest::newRow("water") << "O" << "H2O" << 17.35;
    QTest::newRow("morpholine,n-methyl") << "CN1CCOCC1" << "C5H11NO" << 102.47;
    QTest::newRow("morpholine") << "C1COCCN1" << "C4H9NO" << 85.17;
    QTest::newRow("nitrogen") << "[N][N]" << "N2" << 25.28;
    QTest::newRow("butanethiol") << "CCCCS" << "C4H10S" << 96.25;
    QTest::newRow("dimethylsulfide") << "CSC" << "C2H6S" << 61.66;
    QTest::newRow("dimethyldisulfide") << "CSSC" << "C2H6S2" << 80.17;
    QTest::newRow("carbondisulfide") << "S=C=S" << "CS2" << 57.60;
    QTest::newRow("diethylsulfide") << "CCSCC" << "C4H10S" << 96.25;
    QTest::newRow("n,n-dimethylbenzenesulfonamide") << "CN(C)S(=O)(=O)c1ccccc1" << "C8H11NO2S" << 162.84;
    QTest::newRow("methylsulfone") << "CS(C)(=O)=O" << "C2H6O2S" << 79.24;
    QTest::newRow("diethyldisulfide") << "CCSSCC" << "C4H10S2" << 114.76;
    QTest::newRow("hydrogen sulfide") << "S" << "H2S" << 27.07;
    QTest::newRow("3-methylthiophene") << "Cc1ccsc1" << "C5H6S" << 85.02;
    QTest::newRow("sulfur dioxide") << "O=S=O" << "O2S" << 42.01;
    QTest::newRow("thiophene") << "c1ccsc1" << "C4H4S" << 67.72;
    QTest::newRow("thiophenol") << "Sc1ccccc1" << "C6H6S" << 99.68;
//    QTest::newRow("tetradecane") << "CCCCCCCCCCCCCC" << "C14H30" << 223.90;
    QTest::newRow("1-methylnaphthalene") << "Cc1cccc2ccccc12" << "C11H10" << 139.12;
    QTest::newRow("2-methylnaphthalene") << "Cc2ccc1ccccc1c2" << "C11H10" << 139.12;
    QTest::newRow("1-ethylnaphthalene") << "CCc1cccc2ccccc12" << "C12H12" << 156.41;
    QTest::newRow("2-ethylnaphthalene") << "CCc2ccc1ccccc1c2" << "C12H12" << 156.41;
    QTest::newRow("1,3-dimethylnaphthalene") << "Cc2cc(C)c1ccccc1c2" << "C12H12" << 156.41;
    QTest::newRow("1,4-dimethylnaphthalene") << "Cc1ccc(C)c2ccccc12" << "C12H12" << 156.41;
    QTest::newRow("1,5-dimethylnaphthalene") << "Cc1cccc2c(C)cccc12" << "C12H12" << 156.41;
    QTest::newRow("2,3-dimethylnaphthalene") << "Cc2cc1ccccc1cc2C" << "C12H12" << 156.41;
//    QTest::newRow("acenaphthene") << "C1Cc2cccc3cccc1c23" << "C12H10" << 133.16;
    QTest::newRow("benzo[k]fluoranthene") << "c2ccc1cc3c(cc1c2)c4cccc5cccc3c45" << "C20H12" << 211.83;
    QTest::newRow("chrysene") << "c1ccc3c(c1)ccc4c2ccccc2ccc34" << "C18H12" << 203.13;
    QTest::newRow("triphenylene") << "c1ccc3c(c1)c2ccccc2c4ccccc34" << "C18H12" << 203.13;
    QTest::newRow("acenaphthylene") << "c13cccc2c3c(ccc1)cc2" << "C12H8" << 130.52;
    QTest::newRow("pyrene") << "c1cc2ccc3cccc4ccc(c1)c2c34" << "C16H10" << 171.17;
    QTest::newRow("perylene") << "c1cc2cccc3c4cccc5cccc(c(c1)c23)c45" << "C20H12" << 211.83;
    QTest::newRow("i-hexane") << "CCCC(C)C" << "C6H14" << 112.33;
    QTest::newRow("3-methylpentane") << "CCC(C)CC" << "C6H14" << 112.33;
    QTest::newRow("2,3-dimethylbutane") << "CC(C)C(C)C" << "C6H14" << 112.33;
    QTest::newRow("2-methylhexane") << "CCCCC(C)C" << "C7H16" << 129.63;
    QTest::newRow("2,3-dimethylpentane") << "CCC(C)C(C)C" << "C7H16" << 129.63;
    QTest::newRow("2,4-dimethylpentane") << "CC(C)CC(C)C" << "C7H16" << 129.63;
    QTest::newRow("3-methylheptane") << "CCCCC(C)CC" << "C8H18" << 146.92;
    QTest::newRow("2-methyloctane") << "CC(C)CCCCCC" << "C9H20" << 164.22;
    QTest::newRow("4-methyloctane") << "CCCC(C)CCCC" << "C9H20" << 164.22;
    QTest::newRow("neohexane") << "CCC(C)(C)C" << "C6H14" << 112.33;
    QTest::newRow("3,3-dimethylpentane") << "CCC(C)(C)CC" << "C7H16" << 129.63;
    QTest::newRow("3,3-diethylpentane") << "CCC(CC)(CC)CC" << "C9H20" << 164.22;
    QTest::newRow("2,2,3-trimethylpentane") << "CC(C)(C)C(C)CC" << "C8H18" << 146.92;
    QTest::newRow("isooctane") << "CCCCCC(C)C" << "C8H18" << 146.92;
    QTest::newRow("2,2,5-trimethylhexane") << "CC(C)(C)CCC(C)C" << "C9H20" << 164.22;
    QTest::newRow("cycloheptane") << "C1CCCCCC1" << "C7H14" << 117.27;
    QTest::newRow("methylcyclohexane") << "CC1CCCCC1" << "C7H14" << 117.27;
    QTest::newRow("1,2-dimethylcyclohexane-cis") << "CC1CCCCC1C" << "C8H16" << 134.57;
    QTest::newRow("pentylcyclopentane") << "CCCCCC1CCCC1" << "C10H20" << 169.16;
    QTest::newRow("2-butene-cis") << "CC=CC" << "C4H8" << 75.10;
    QTest::newRow("1-heptene") << "CCCCCC=C" << "C7H14" << 126.99;
    QTest::newRow("2-heptene-tr") << "CCCCC=CC" << "C7H14" << 126.99;
    QTest::newRow("1-nonene") << "CCCCCCCC=C" << "C9H18" << 161.58;
    QTest::newRow("heptane") << "CCCCCCC" << "C7H16" << 129.63;
    QTest::newRow("1,5-hexadiene") << "C=CCCC=C" << "C6H10" << 107.06;
    QTest::newRow("cyclopentadiene") << "C1C=CC=C1" << "C5H6" << 77.41;
    QTest::newRow("myrcene") << "C=CC(=C)CCC=C(C)C" << "C10H16" << 173.61;
    QTest::newRow("1-heptyne") << "CCCCCC#C" << "C7H12" << 124.36;
    QTest::newRow("1-nonyne") << "CCCCCCCC#C" << "C9H16" << 158.95;
    QTest::newRow("nonane") << "CCCCCCCCC" << "C9H20" << 164.22;
    QTest::newRow("propylbenzene") << "CCCc1ccccc1" << "C9H12" << 133.05;
    QTest::newRow("1-ethyl-2-methylbenzene") << "CCc1ccccc1C" << "C9H12" << 133.05;
    QTest::newRow("benzene,1-ethyl-4-methyl") << "CCc1ccc(C)cc1" << "C9H12" << 133.05;
    QTest::newRow("1,2,3-trimethylbenzene") << "Cc1cccc(C)c1C" << "C9H12" << 133.05;
    QTest::newRow("1,2,4-trimethylbenzene") << "Cc1ccc(C)c(C)c1" << "C9H12" << 133.05;
    QTest::newRow("1,3,5-trimethylbenzene") << "Cc1cc(C)cc(C)c1" << "C9H12" << 133.05;
    QTest::newRow("decane") << "CCCCCCCCCC" << "C10H22" << 181.52;
    QTest::newRow("butylbenzene") << "CCCCc1ccccc1" << "C10H14" << 150.35;
    QTest::newRow("1,2-diethylbenzene") << "CCc1ccccc1CC" << "C10H14" << 150.35;
    QTest::newRow("benzene,1,3-diethyl") << "CCc1cccc(CC)c1" << "C10H14" << 150.35;
    QTest::newRow("1,4-diethylbenzene") << "CCc1ccc(CC)cc1" << "C10H14" << 150.35;
    QTest::newRow("pentylbenzene") << "CCCCCc1ccccc1" << "C11H16" << 167.65;
    QTest::newRow("pentamethylbenzene") << "Cc1cc(C)c(C)c(C)c1C" << "C11H16" << 167.65;
    QTest::newRow("hexamethylbenzene") << "Cc1c(C)c(C)c(C)c(C)c1C" << "C12H18" << 184.94;
    QTest::newRow("heptylbenzene") << "CCCCCCCc1ccccc1" << "C13H20" << 202.24;
    QTest::newRow("octylbenzene") << "CCCCCCCCc1ccccc1" << "C14H22" << 219.53;
    QTest::newRow("dodecane") << "CCCCCCCCCCCC" << "C12H26" << 216.11;
    QTest::newRow("benzene,1-methyl-2-i-propyl") << "CC(C)c1ccccc1C" << "C10H14" << 150.35;
    QTest::newRow("p-i-propyltoluene") << "CC(C)c1ccc(C)cc1" << "C10H14" << 150.35;
    QTest::newRow("indane") << "C2Cc1ccccc1C2" << "C9H10" << 120.70;
    QTest::newRow("3-methylstyrene") << "c1ccc(C)cc1C=C" << "C9H10" << 130.42;
    QTest::newRow("4-methylstyrene") << "c1cc(C)ccc1C=C" << "C9H10" << 130.42;
    QTest::newRow("fluorene") << "C2c1ccccc1c3ccccc23" << "C13H10" << 158.72;
    QTest::newRow("n-octadecan-1-ol") << "CCCCCCCCCCCCCCCCCCO" << "C18H38O" << 328.67;
    QTest::newRow("2-butenal-tr") << "CC=CC=O" << "C4H6O" << 81.26;
    QTest::newRow("2-hexenal-tr") << "CCCC=CC=O" << "C6H10O" << 115.85;
    QTest::newRow("2-methylpropenal") << "CC(=C)C=O" << "C4H6O" << 81.26;
    QTest::newRow("2-ethyl-2-hexenal") << "CCCC=C(CC)C=O" << "C8H14O" << 150.44;
    QTest::newRow("2-hexanone") << "CCCCC(C)=O" << "C6H12O" << 118.49;
    QTest::newRow("3-hexanone") << "CCCC(=O)CC" << "C6H12O" << 118.49;
    QTest::newRow("3-heptanone") << "CCCCC(=O)CC" << "C7H14O" << 135.78;
    QTest::newRow("octan-3-one") << "CCCCCC(=O)CC" << "C8H16O" << 153.08;
//    QTest::newRow("5-nonanone") << "CCCCC(=O)CCCC" << "C9H18O" << 155.71;
    QTest::newRow("2-nonanone") << "CCCCCCCC(C)=O" << "C9H18O" << 170.37;
    QTest::newRow("2-decanone") << "CCCCCCCCC(C)=O" << "C10H20O" << 187.67;
    QTest::newRow("2-methyl-pentan-1-ol") << "CC(C)CCCO" << "C6H14O" << 121.12;
    QTest::newRow("4-methyl-2-pentanone") << "CC(C)CC(C)=O" << "C6H12O" << 118.49;
    QTest::newRow("2,4-dimethyl-3-pentanone") << "CC(C)C(=O)C(C)C" << "C7H14O" << 135.78;
    QTest::newRow("5-methyl-3-heptanone") << "CCC(=O)CC(C)CC" << "C8H16O" << 153.08;
    QTest::newRow("2,6-dimethyl-4-heptanone") << "CC(C)CC(=O)CC(C)C" << "C9H18O" << 170.37;
    QTest::newRow("2-methyl-cyclohexanone") << "C1CCCC(C)C1=O" << "C7H12O" << 123.43;
//    QTest::newRow("1-buten-3-one") << "CCC(=O)C" << "C4H8O" << 81.26;
    QTest::newRow("p-methylacetophenone") << "CC(=O)c1ccc(C)cc1" << "C9H10O" << 139.21;
//    QTest::newRow("anthraquinone") << "O=C2c1ccccc1C(=O)c3ccccc23" << "C14H8O2" << 188.32;
    QTest::newRow("formicacid") << "OC=O" << "CH2O2" << 40.80;
    QTest::newRow("1-propanol,2,2-dimethyl") << "CC(C)(C)CO" << "C5H12O" << 103.83;
    QTest::newRow("hexanoicacid") << "CCCCCC(O)=O" << "C6H12O2" << 127.28;
    QTest::newRow("heptanoicacid") << "CCCCCCC(O)=O" << "C7H14O2" << 144.57;
    QTest::newRow("octanoicacid") << "CCCCCCCC(O)=O" << "C8H16O2" << 161.87;
    QTest::newRow("i-valericacid") << "CC(C)CC(O)=O" << "C5H10O2" << 109.98;
    QTest::newRow("2-methylbutyric acid") << "CCC(C)C(O)=O" << "C5H10O2" << 109.98;
    QTest::newRow("2-ethyl-butanoicacid") << "CCC(CC)C(O)=O" << "C6H12O2" << 127.28;
    QTest::newRow("formicacid,butylester") << "CCCCOC=O" << "C5H10O2" << 109.98;
    QTest::newRow("i-butylformate") << "CC(C)COC=O" << "C5H10O2" << 109.98;
    QTest::newRow("2-pentanol") << "CCCC(C)O" << "C5H12O" << 103.83;
    QTest::newRow("sec-butylacetate") << "CC(=O)OC(C)CC" << "C6H12O2" << 127.28;
//    QTest::newRow("methyl-tert-butylacetate") << "CC(C)(C)(=O)OC" << "C5H12O2" << 127.28;
    QTest::newRow("hexanoicacid,methylester") << "CCCCCC(=O)OC" << "C7H14O2" << 144.57;
    QTest::newRow("butylpropionate") << "CCCCOC(=O)CC" << "C7H14O2" << 144.57;
    QTest::newRow("aceticacid,hexylester") << "CCCCCCOC(C)=O" << "C8H16O2" << 161.87;
    QTest::newRow("methyloctanoate") << "CCCCCCCC(=O)OC" << "C9H18O2" << 179.16;
    QTest::newRow("tert-butylacetate") << "CC(=O)OC(C)(C)C" << "C6H12O2" << 127.28;
    QTest::newRow("4-methylpentylacetate") << "CC(=O)OCCC(C)CC" << "C8H16O2" << 161.87;
    QTest::newRow("heptan-3-ol") << "CCCCC(O)CC" << "C7H16O" << 138.42;
    QTest::newRow("propiolactone") << "O=C1CCO1" << "C3H4O2" << 63.03;
    QTest::newRow("vinylacetate") << "CC(=O)OC=C" << "C4H6O2" << 90.05;
    QTest::newRow("methacrylicacid,methylester") << "COC(=O)C(C)=C" << "C5H8O2" << 107.34;
    QTest::newRow("methacrylicacid,i-butylester") << "CC(C)COC(=O)C(C)=C" << "C8H14O2" << 159.23;
    QTest::newRow("oxalic acid,dimethyl ester") << "COC(=O)C(=O)OC" << "C4H6O4" << 107.63;
    QTest::newRow("2-nonanol") << "CCCCCCCC(O)C" << "C9H20O" << 173.01;
    QTest::newRow("butanedioicacid,methylester") << "COC(=O)CCC(=O)OC" << "C6H10O4" << 142.22;
    QTest::newRow("butanedioicacid,ethylester") << "CCOC(=O)CCC(=O)OCC" << "C8H14O4" << 176.81;
    QTest::newRow("ethylbenzoate") << "CCOC(=O)c1ccccc1" << "C9H10O2" << 148.00;
    QTest::newRow("diethylphthalate") << "CCOC(=O)c1ccccc1C(=O)OCC" << "C12H14O4" << 214.83;
    QTest::newRow("dipropylphthalate") << "CCCOC(=O)c1ccccc1C(=O)OCCC" << "C14H18O4" << 249.42;
    QTest::newRow("dihexylphthalate") << "CCCCCCOC(=O)c1ccccc1C(=O)OCCCCCC" << "C20H30O4" << 353.20;
    QTest::newRow("3-methyl-2-butanol") << "CC(C)C(C)O" << "C5H12O" << 103.83;
    QTest::newRow("propylenecarbonate") << "CC1COC(=O)O1" << "C4H6O3" << 89.12;
    QTest::newRow("aceticanhydride") << "CC(=O)OC(=O)C" << "C4H6O3" << 98.84;
//    QTest::newRow("maleicanhydride") << "C1(=O)C=CC(=O)O1" << "C4H2O3" << 83.85;
    QTest::newRow("i-propoxyethanol") << "CC(C)OCCO" << "C5H12O2" << 112.62;
    QTest::newRow("2-phenoxyethanol") << "OCCOc1ccccc1" << "C8H10O2" << 133.34;
//    QTest::newRow("4-methyl-2-pentanol") << "CC(C)CC(C)O" << "C6H14O" << 107.94;
    QTest::newRow("m-methoxyphenol") << "COc1cccc(O)c1" << "C7H8O2" << 116.04;
    QTest::newRow("2,6-dimethoxyphenol") << "COc1cccc(OC)c1O" << "C8H10O3" << 142.13;
    QTest::newRow("o-hydroxybenzaldehyde") << "Oc1ccccc1C=O" << "C7H6O2" << 113.41;
    QTest::newRow("p-methoxyacetophenone") << "COc1ccc(cc1)C(C)=O" << "C9H10O2" << 148.00;
    QTest::newRow("cyclopentanol") << "OC1CCCC1" << "C5H10O" << 91.47;
    QTest::newRow("ethylacetoacetate") << "CCOC(=O)CC(C)=O" << "C6H10O3" << 133.43;
    QTest::newRow("cycloheptanol") << "OC1CCCCCC1" << "C7H14O" << 126.06;
    QTest::newRow("cyclohexanol,2-methyl") << "C1CCCC(C)C1O" << "C7H14O" << 126.06;
    QTest::newRow("cyclohexanol,3-methyl") << "C1CCC(C)CC1O" << "C7H14O" << 126.06;
    QTest::newRow("cyclohexanol,4-methyl") << "CC1CCC(O)CC1" << "C7H14O" << 126.06;
    QTest::newRow("a-terpineol") << "CC1=CCC(CC1)C(C)(C)O" << "C10H18O" << 175.31;
    QTest::newRow("phenylethyl alcohol") << "OCCc1ccccc1" << "C8H10O" << 124.55;
    QTest::newRow("2,3-dimethylphenol") << "Cc1cccc(O)c1C" << "C8H10O" << 124.55;
    QTest::newRow("2,5-dimethylphenol") << "Cc1ccc(C)c(O)c1" << "C8H10O" << 124.55;
    QTest::newRow("2,6-dimethylphenol") << "Cc1cccc(C)c1O" << "C8H10O" << 124.55;
    QTest::newRow("3,4-dimethylphenol") << "Cc1ccc(O)cc1C" << "C8H10O" << 124.55;
    QTest::newRow("3,5-dimethylphenol") << "Cc1cc(C)cc(O)c1" << "C8H10O" << 124.55;
    QTest::newRow("p-ethylphenol") << "CCc1ccc(O)cc1" << "C8H10O" << 124.55;
    QTest::newRow("p-propylphenol") << "CCCc1ccc(O)cc1" << "C9H12O" << 141.84;
    QTest::newRow("1-naphthol") << "Oc1cccc2ccccc12" << "C10H8O" << 130.61;
    QTest::newRow("2-naphthol") << "Oc2ccc1ccccc1c2" << "C10H8O" << 130.61;
    QTest::newRow("m-dihydroxybenzene") << "Oc1cccc(O)c1" << "C6H6O2" << 98.75;
    QTest::newRow("methoxyisopropane") << "CC(OC)C" << "C4H10O" << 86.53;
    QTest::newRow("methoxybutane") << "CCCCOC" << "C5H12O" << 103.83;
    QTest::newRow("butoxyvinyl") << "CCCCOC=C" << "C6H12O" << 118.49;
    QTest::newRow("iso-butoxyvinyl") << "C(C)CCOC=C" << "C6H12O" << 118.49;
    QTest::newRow("vinyl ether") << "C=COC=C" << "C4H6O" << 81.26;
    QTest::newRow("1,2-dibutoxyethane") << "CCCCOCCOCCCC" << "C10H22O2" << 199.10;
    QTest::newRow("eucalyptol") << "CC12CCC(CC1)C(C)(C)O2" << "C10H18O" << 165.59;
    QTest::newRow("2-methylfuran") << "Cc1ccco1" << "C5H6O" << 75.30;
    QTest::newRow("dibenzofuran") << "c1ccc2c(c1)oc3ccccc23" << "C12H8O" << 139.31;
    QTest::newRow("1-dodecanol") << "CCCCCCCCCCCCO" << "C12H26O" << 224.90;
    QTest::newRow("heptanal") << "CCCCCCC=O" << "C7H14O" << 135.78;
    QTest::newRow("nonanal") << "CCCCCCCCC=O" << "C9H18O" << 170.37;
    QTest::newRow("2-ethylbutanal") << "CCC(CC)C=O" << "C6H12O" << 118.49;
    QTest::newRow("proargylalcohol") << "C#CCO" << "C3H4O" << 63.96;
    QTest::newRow("2-methyl-1,3-butadiene") << "CC(=C)C=C" << "C5H8" << 89.76;
    QTest::newRow("methylformamide") << "CNC=O" << "C2H5NO" << 60.30;
    QTest::newRow("anisole") << "COc1ccccc1" << "C7H8O" << 107.25;
//    QTest::newRow("benzoicacid,methylester") << "COC(=O)c1ccccc1" << "C8H8O2" << 145.36;
    QTest::newRow("3-methylindole") << "Cc1c[nH]c2ccccc12" << "C9H9N" << 118.16;
    QTest::newRow("diethylcarbonate") << "CCOC(=O)OCC" << "C5H10O3" << 118.77;
    QTest::newRow("dimethylcarbonate") << "COC(=O)OC" << "C3H6O3" << 84.18;
    QTest::newRow("aceticacid,methylester") << "COC(C)=O" << "C3H6O2" << 75.39;
    QTest::newRow("benzene,trifluoromethyl") << "FC(F)(F)c1ccccc1" << "C7H5F3" << 116.67;
    QTest::newRow("dibenzothiophene") << "c1ccc2c(c1)sc3ccccc23" << "C12H8S" << 149.03;
//    QTest::newRow("methylthiobenzene") << "CSc1ccccc1" << "C7H8S" << 122.70;
    QTest::newRow("benzene,ethynyl") << "C#Cc1ccccc1" << "C8H6" << 110.49;
    QTest::newRow("isophorone") << "C1C(C)(C)CC(C)=CC1=O" << "C9H14O" << 155.38;
    QTest::newRow("phenol,2-methoxy-4-allyl") << "COc1cc(CC=C)ccc1O" << "C10H12O2" << 165.29;
    QTest::newRow("salicylic acid") << "OC(=O)c1ccccc1O" << "C7H6O3" << 122.20;
    QTest::newRow("furfural") << "O=Cc1ccco1" << "C5H4O2" << 81.45;
    QTest::newRow("diphenylether") << "O(c1ccccc1)c2ccccc2" << "C12H10O" << 162.57;
    QTest::newRow("paraldehyde") << "CC1OC(C)OC(C)O1" << "C6H12O3" << 126.35;
//    QTest::newRow("34dinitrophenol") << "Oc1ccc(N(=O)=O)c(c1)N(=O)=O" << "C6H4N2O5" << 162.46;
    QTest::newRow("resorcinol") << "Oc1cccc(O)c1" << "C6H6O2" << 98.75;
    QTest::newRow("hydroquinone") << "Oc1ccc(O)cc1" << "C6H6O2" << 98.75;
    QTest::newRow("methylgallate") << "COC(=O)c1cc(O)c(O)c(O)c1" << "C8H8O5" << 157.07;
    QTest::newRow("acetaminophen") << "CC(=O)Nc1ccc(O)cc1" << "C8H9NO2" << 141.70;
    QTest::newRow("123trihydroxybenzene") << "Oc1cccc(O)c1O" << "C6H6O3" << 107.54;
    QTest::newRow("135trihydroxybenzene") << "Oc1cc(O)cc(O)c1" << "C6H6O3" << 107.54;
    QTest::newRow("ophthalicacid") << "OC(=O)c1ccccc1C(=O)O" << "C8H6O4" << 145.65;
    QTest::newRow("phydroxybenzaldehyde") << "Oc1ccc(C=O)cc1" << "C7H6O2" << 113.41;
    QTest::newRow("mcyanophenol") << "Oc1cccc(C#N)c1" << "C7H5NO" << 112.98;
    QTest::newRow("phenytoin") << "O=C1NC(=O)C(N1)(c2ccccc2)c3ccccc3" << "C15H12N2O2" << 227.61;
    QTest::newRow("tetraethylurea") << "CCN(CC)C(=O)N(CC)CC" << "C9H20N2O" << 192.37;
    QTest::newRow("triethylphosphate") << "CCOP(=O)(OCC)OCC" << "C6H15O4P" << 167.32;
//    QTest::newRow("COP(O)(NC(C)O)SC") << "COP(O)(NC(C)O)SC" << "C4H14NO3PS" << 156.08;
    QTest::newRow("ephedrine") << "c1ccccc1C(O)C(C)NC" << "C10H15NO" << 170.14;
    QTest::newRow("triphenylphosphineoxide") << "O=P(c1ccccc1)(c2ccccc2)c3ccccc3" << "C18H15OP" << 255.00;
    QTest::newRow("nicotine") << "CN1CCCC1c2cccnc2" << "C10H14N2" << 159.99;
    QTest::newRow("meperidine") << "CCOC(=O)C1(CCN(C)CC1)c2ccccc2" << "C15H21NO2" << 250.41;
    QTest::newRow("antipyrine") << "Cc1cc(=O)n(c2ccccc2)n1C" << "C11H12N2O" << 169.90;
    QTest::newRow("3bromoacetanilide") << "CC(=O)Nc1cccc(Br)c1" << "C8H8BrNO" << 152.19;
    QTest::newRow("24dinitrophenol") << "Oc1ccc(cc1N(=O)=O)N(=O)=O" << "C6H4N2O5" << 141.84;
    QTest::newRow("gallicacid") << "OC(=O)c1cc(O)c(O)c(O)c1" << "C7H6O5" << 139.78;
    QTest::newRow("succinicacid") << "OC(=O)CCC(=O)O" << "C4H6O4" << 107.63;
    QTest::newRow("glyceryltriacetate") << "CC(=O)OCC(COC(=O)C)OC(=O)C" << "C9H14O6" << 209.05;
    QTest::newRow("1methyl2cyanoguanidine") << "CNC(=NC#N)N" << "C3H6N4" << 96.52;
    QTest::newRow("dimethylsulfate") << "COS(=O)(=O)OC" << "C2H6O4S" << 96.82;
    QTest::newRow("hexamethylphosphotriamide") << "CN(C)P(N(C)C)N(C)C" << "C6H18N3P" << 165.15;
    QTest::newRow("methadone") << "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c2ccccc2" << "C21H27NO" << 326.59;
    QTest::newRow("nitrofurazone") << "NC(=O)NN=Cc1ccc(o1)N(=O)=O" << "C6H6N4O4" << 155.04;
//    QTest::newRow("thymine") << "CC1=CNC(=O)NC1=O" << "C5H6N2O2" << 114.34;
    QTest::newRow("barbital") << "O=C1NC(=O)NC(=O)C1(CC)CC" << "C8H12N2O3" << 175.02;
    QTest::newRow("pentobarbital") << "O=C1NC(=O)NC(=O)C1(CC)C(C)CCC" << "C11H18N2O3" << 226.91;
}

void VabcTest::vabc()
{
    QFETCH(QString, smilesString);
    QFETCH(QString, formulaString);
    QFETCH(double, vabc);

    QByteArray smiles = smilesString.toAscii();
    QByteArray formula = formulaString.toAscii();

    chemkit::Molecule molecule(smiles.constData(), "smiles");
    QCOMPARE(molecule.formula().c_str(), formula.constData());
    QCOMPARE(qRound(molecule.descriptor("vabc").toDouble()), qRound(vabc));
}

QTEST_APPLESS_MAIN(VabcTest)
