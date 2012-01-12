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

#include "pubchemfingerprint.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

// PubChem Fingerprint Specification:
// ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt

PubChemFingerprint::PubChemFingerprint()
    : chemkit::Fingerprint("pubchem")
{
}

PubChemFingerprint::~PubChemFingerprint()
{
}

chemkit::Bitset PubChemFingerprint::value(const chemkit::Molecule *molecule) const
{
    chemkit::Bitset bitset(881);

    // section 1 - hierarchic element counts
    size_t hydrogenCount = molecule->atomCount(chemkit::Atom::Hydrogen);
    bitset[0] = hydrogenCount >= 4;
    bitset[1] = hydrogenCount >= 8;
    bitset[2] = hydrogenCount >= 16;
    bitset[3] = hydrogenCount >= 32;

    size_t lithiumCount = molecule->atomCount(chemkit::Atom::Lithium);
    bitset[4] = lithiumCount >= 1;
    bitset[5] = lithiumCount >= 2;

    size_t boronCount = molecule->atomCount(chemkit::Atom::Boron);
    bitset[6] = boronCount >= 1;
    bitset[7] = boronCount >= 2;
    bitset[8] = boronCount >= 4;

    size_t carbonCount = molecule->atomCount(chemkit::Atom::Carbon);
    bitset[9] = carbonCount >= 2;
    bitset[10] = carbonCount >= 4;
    bitset[11] = carbonCount >= 8;
    bitset[12] = carbonCount >= 16;
    bitset[13] = carbonCount >= 32;

    size_t nitrogenCount = molecule->atomCount(chemkit::Atom::Nitrogen);
    bitset[14] = nitrogenCount >= 1;
    bitset[15] = nitrogenCount >= 2;
    bitset[16] = nitrogenCount >= 4;
    bitset[17] = nitrogenCount >= 8;

    size_t oxygenCount = molecule->atomCount(chemkit::Atom::Oxygen);
    bitset[18] = oxygenCount >= 1;
    bitset[19] = oxygenCount >= 2;
    bitset[20] = oxygenCount >= 4;
    bitset[21] = oxygenCount >= 8;
    bitset[22] = oxygenCount >= 16;

    size_t fluorineCount = molecule->atomCount(chemkit::Atom::Fluorine);
    bitset[23] = fluorineCount >= 1;
    bitset[24] = fluorineCount >= 2;
    bitset[25] = fluorineCount >= 4;

    size_t sodiumCount = molecule->atomCount(chemkit::Atom::Sodium);
    bitset[26] = sodiumCount >= 1;
    bitset[27] = sodiumCount >= 2;

    size_t siliconCount = molecule->atomCount(chemkit::Atom::Silicon);
    bitset[28] = siliconCount >= 1;
    bitset[29] = siliconCount >= 2;

    size_t phosphorusCount = molecule->atomCount(chemkit::Atom::Phosphorus);
    bitset[30] = phosphorusCount >= 1;
    bitset[31] = phosphorusCount >= 2;
    bitset[32] = phosphorusCount >= 4;

    size_t sulfurCount = molecule->atomCount(chemkit::Atom::Sulfur);
    bitset[33] = sulfurCount >= 1;
    bitset[34] = sulfurCount >= 2;
    bitset[35] = sulfurCount >= 4;
    bitset[36] = sulfurCount >= 8;

    size_t chlorineCount = molecule->atomCount(chemkit::Atom::Chlorine);
    bitset[37] = chlorineCount >= 1;
    bitset[38] = chlorineCount >= 2;
    bitset[39] = chlorineCount >= 4;
    bitset[40] = chlorineCount >= 8;

    size_t potassiumCount = molecule->atomCount(chemkit::Atom::Potassium);
    bitset[41] = potassiumCount >= 1;
    bitset[42] = potassiumCount >= 2;

    size_t bromineCount = molecule->atomCount(chemkit::Atom::Bromine);
    bitset[43] = bromineCount >= 1;
    bitset[44] = bromineCount >= 2;
    bitset[45] = bromineCount >= 4;

    size_t iodineCount = molecule->atomCount(chemkit::Atom::Iodine);
    bitset[46] = iodineCount >= 1;
    bitset[47] = iodineCount >= 2;
    bitset[48] = iodineCount >= 4;

    bitset[49] = molecule->contains(chemkit::Atom::Beryllium);
    bitset[50] = molecule->contains(chemkit::Atom::Magnesium);
    bitset[51] = molecule->contains(chemkit::Atom::Aluminum);
    bitset[52] = molecule->contains(chemkit::Atom::Calcium);
    bitset[53] = molecule->contains(chemkit::Atom::Scandium);
    bitset[54] = molecule->contains(chemkit::Atom::Titanium);
    bitset[55] = molecule->contains(chemkit::Atom::Vanadium);
    bitset[56] = molecule->contains(chemkit::Atom::Chromium);
    bitset[57] = molecule->contains(chemkit::Atom::Manganese);
    bitset[58] = molecule->contains(chemkit::Atom::Iron);
    bitset[59] = molecule->contains(chemkit::Atom::Cobalt);
    bitset[60] = molecule->contains(chemkit::Atom::Nickel);
    bitset[61] = molecule->contains(chemkit::Atom::Copper);
    bitset[62] = molecule->contains(chemkit::Atom::Zinc);
    bitset[63] = molecule->contains(chemkit::Atom::Gallium);
    bitset[64] = molecule->contains(chemkit::Atom::Germanium);
    bitset[65] = molecule->contains(chemkit::Atom::Arsenic);
    bitset[66] = molecule->contains(chemkit::Atom::Selenium);
    bitset[67] = molecule->contains(chemkit::Atom::Krypton);
    bitset[68] = molecule->contains(chemkit::Atom::Rubidium);
    bitset[69] = molecule->contains(chemkit::Atom::Strontium);
    bitset[70] = molecule->contains(chemkit::Atom::Yttrium);
    bitset[71] = molecule->contains(chemkit::Atom::Zirconium);
    bitset[72] = molecule->contains(chemkit::Atom::Niobium);
    bitset[73] = molecule->contains(chemkit::Atom::Molybdenum);
    bitset[74] = molecule->contains(chemkit::Atom::Ruthenium);
    bitset[75] = molecule->contains(chemkit::Atom::Rhodium);
    bitset[76] = molecule->contains(chemkit::Atom::Palladium);
    bitset[77] = molecule->contains(chemkit::Atom::Silver);
    bitset[78] = molecule->contains(chemkit::Atom::Cadmium);
    bitset[79] = molecule->contains(chemkit::Atom::Indium);
    bitset[80] = molecule->contains(chemkit::Atom::Tin);
    bitset[81] = molecule->contains(chemkit::Atom::Antimony);
    bitset[82] = molecule->contains(chemkit::Atom::Tellurium);
    bitset[83] = molecule->contains(chemkit::Atom::Xenon);
    bitset[84] = molecule->contains(chemkit::Atom::Cesium);
    bitset[85] = molecule->contains(chemkit::Atom::Barium);
    bitset[86] = molecule->contains(chemkit::Atom::Lutetium);
    bitset[87] = molecule->contains(chemkit::Atom::Hafnium);
    bitset[88] = molecule->contains(chemkit::Atom::Tantalum);
    bitset[89] = molecule->contains(chemkit::Atom::Tungsten);
    bitset[90] = molecule->contains(chemkit::Atom::Rhenium);
    bitset[91] = molecule->contains(chemkit::Atom::Osmium);
    bitset[92] = molecule->contains(chemkit::Atom::Iridium);
    bitset[93] = molecule->contains(chemkit::Atom::Platinum);
    bitset[94] = molecule->contains(chemkit::Atom::Gold);
    bitset[95] = molecule->contains(chemkit::Atom::Mercury);
    bitset[96] = molecule->contains(chemkit::Atom::Thallium);
    bitset[97] = molecule->contains(chemkit::Atom::Lead);
    bitset[98] = molecule->contains(chemkit::Atom::Bismuth);
    bitset[99] = molecule->contains(chemkit::Atom::Lanthanum);
    bitset[100] = molecule->contains(chemkit::Atom::Cerium);
    bitset[101] = molecule->contains(chemkit::Atom::Praseodymium);
    bitset[102] = molecule->contains(chemkit::Atom::Neodymium);
    bitset[103] = molecule->contains(chemkit::Atom::Promethium);
    bitset[104] = molecule->contains(chemkit::Atom::Samarium);
    bitset[105] = molecule->contains(chemkit::Atom::Europium);
    bitset[106] = molecule->contains(chemkit::Atom::Gadolinium);
    bitset[107] = molecule->contains(chemkit::Atom::Terbium);
    bitset[108] = molecule->contains(chemkit::Atom::Dysprosium);
    bitset[109] = molecule->contains(chemkit::Atom::Holmium);
    bitset[110] = molecule->contains(chemkit::Atom::Erbium);
    bitset[111] = molecule->contains(chemkit::Atom::Thulium);
    bitset[112] = molecule->contains(chemkit::Atom::Ytterbium);
    bitset[113] = molecule->contains(chemkit::Atom::Technetium);
    bitset[114] = molecule->contains(chemkit::Atom::Uranium);

    // section 2 - ring counts
    // TODO

    // section 3 - simple atom pairs
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Hydrogen)){
            bitset[263] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Lithium)){
            bitset[264] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Boron)){
            bitset[265] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Carbon)){
            bitset[266] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Oxygen)){
            bitset[267] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Fluorine)){
            bitset[268] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Phosphorus)){
            bitset[269] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Sulfur)){
            bitset[270] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Lithium, chemkit::Atom::Chlorine)){
            bitset[271] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Hydrogen)){
            bitset[272] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Boron)){
            bitset[273] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Carbon)){
            bitset[274] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Nitrogen)){
            bitset[275] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Oxygen)){
            bitset[276] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Fluorine)){
            bitset[277] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Silicon)){
            bitset[278] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Phosphorus)){
            bitset[279] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Sulfur)){
            bitset[280] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Chlorine)){
            bitset[281] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Boron, chemkit::Atom::Bromine)){
            bitset[282] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Hydrogen)){
            bitset[283] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Carbon)){
            bitset[284] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Nitrogen)){
            bitset[285] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Oxygen)){
            bitset[286] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Fluorine)){
            bitset[287] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Sodium)){
            bitset[288] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Magnesium)){
            bitset[289] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Aluminum)){
            bitset[290] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Silicon)){
            bitset[291] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Phosphorus)){
            bitset[292] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Sulfur)){
            bitset[293] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Chlorine)){
            bitset[294] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Arsenic)){
            bitset[295] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Selenium)){
            bitset[296] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Bromine)){
            bitset[297] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Carbon, chemkit::Atom::Iodine)){
            bitset[298] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Hydrogen)){
            bitset[299] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Nitrogen)){
            bitset[300] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Oxygen)){
            bitset[301] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Fluorine)){
            bitset[302] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Silicon)){
            bitset[303] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Phosphorus)){
            bitset[304] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Sulfur)){
            bitset[305] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Chlorine)){
            bitset[306] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Nitrogen, chemkit::Atom::Bromine)){
            bitset[307] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Hydrogen)){
            bitset[308] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Oxygen)){
            bitset[309] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Magnesium)){
            bitset[310] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Sodium)){
            bitset[311] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Aluminum)){
            bitset[312] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Silicon)){
            bitset[313] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Phosphorus)){
            bitset[314] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Oxygen, chemkit::Atom::Potassium)){
            bitset[315] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Fluorine, chemkit::Atom::Phosphorus)){
            bitset[316] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Fluorine, chemkit::Atom::Sulfur)){
            bitset[317] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Aluminum, chemkit::Atom::Hydrogen)){
            bitset[318] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Aluminum, chemkit::Atom::Chlorine)){
            bitset[319] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Silicon, chemkit::Atom::Hydrogen)){
            bitset[320] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Silicon, chemkit::Atom::Silicon)){
            bitset[321] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Silicon, chemkit::Atom::Chlorine)){
            bitset[322] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Phosphorus, chemkit::Atom::Hydrogen)){
            bitset[323] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Phosphorus, chemkit::Atom::Phosphorus)){
            bitset[324] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Arsenic, chemkit::Atom::Hydrogen)){
            bitset[325] = 1;
        }
        else if(bond->containsBoth(chemkit::Atom::Arsenic, chemkit::Atom::Arsenic)){
            bitset[326] = 1;
        }
    }

    // section 4 - simple atom nearest neighbors
    // TODO

    // section 5 - detailed atom neighborhoods
    // TODO

    // section 6 - simples SMARTS patterns
    // TODO

    // section 7 - complex SMARTS patterns
    // TODO

    return bitset;
}
