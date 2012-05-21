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

#include "apoldescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

namespace {

// The following atomic polarizabilities were taken from:
// http://www.sunysccc.edu/academic/mst/ptable/p-table2.htm
const chemkit::Real atomicPolarizabilities[] = {
    0.0,
    0.666793, // Hydrogen
    0.204956, // Helium
    24.3, // Lithium
    5.60, // Beryllium
    3.03, // Boron
    1.76, // Carbon
    1.10, // Nitrogen
    0.802, // Oxygen
    0.557, // Fluorine
    0.3956, // Neon
    24.08, // Sodium
    10.6, // Magnesium
    6.8, // Aluminum
    5.38, // Silicon
    3.63, // Phosphorus
    2.90, // Sulfur
    2.18, // Chlorine
    1.6411, // Argon
    43.4, // Potassium
    22.8, // Calcium
    17.8, // Scandium
    14.6, // Titanium
    12.4, // Vanadium
    11.6, // Chromium
    9.4, // Manganese
    8.4, // Iron
    7.5, // Cobalt
    6.8, // Nickel
    6.1, // Copper
    7.1, // Zinc
    8.12, // Gallium
    6.07, // Germanium
    4.31, // Arsenic
    3.77, // Selenium
    3.05, // Bromine
    2.4844, // Krypton
    47.3, // Rubidium
    27.6, // Strontium
    22.7, // Yttrium
    17.9, // Zirconium
    15.7, // Niobium
    12.8, // Molybdenum
    11.4, // Technetium
    9.6, // Ruthenium
    8.6, // Rhodium
    4.8, // Palladium
    7.2, // Silver
    7.2, // Cadmium
    9.1, // Indium
    7.7, // Tin
    6.6, // Antimony
    5.5, // Tellurium
    4.7, // Iodine
    4.044, // Xenon
    59.6, // Cesium
    39.7, // Barium
    31.1, // Lanthanum
    29.6, // Cerium
    28.2, // Praseodymium
    31.4, // Neodymium
    30.1, // Promethium
    28.8, // Samarium
    22.7, // Europium
    23.5, // Gadolinium
    25.5, // Terbium
    24.5, // Dysprosium
    23.6, // Holmium
    22.7, // Erbium
    21.8, // Thulium
    21.0, // Ytterbium
    21.9, // Lutetium
    16.2, // Hafnium
    13.1, // Tantalum
    11.1, // Tungsten
    9.7, // Rhenium
    8.5, // Osmium
    7.6, // Iridium
    6.5, // Platinum
    5.8, // Gold
    5.7, // Mercury
    7.5, // Thallium
    6.8, // Lead
    7.4, // Bismuth
    6.8, // Polonium
    6.0, // Astatine
    5.3, // Radon
    48.7, // Francium
    38.3, // Radium
    32.1, // Actinium
    32.1, // Thorium
    25.4, // Protactinium
    24.9, // Uranium
    24.8, // Neptunium
    24.5, // Plutonium
    23.3, // Americium
    23.0, // Curium
    22.7, // Berkelium
    20.5, // Californium
    19.7, // Einsteinium
    23.8, // Fermium
    18.2, // Mendelevium
    17.5, // Nobelium
    0, // Lawrencium
    0, // Rutherfordium
    0, // Dubnium
    0, // Seaborgium
    0, // Bohrium
    0, // Hassium
    0 // Meitnerium
};

} // end anonymous namespace

// === ApolDescriptor ====================================================== //
ApolDescriptor::ApolDescriptor()
    : chemkit::MolecularDescriptor("apol")
{
    setDimensionality(1);
}

// Returns the sum of atomic polarizabilities for each atom in the molecule.
chemkit::Variant ApolDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real value = 0;

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        value += atomicPolarizabilities[atom->atomicNumber()];
    }

    return value;
}

// === BpolDescriptor ====================================================== //
BpolDescriptor::BpolDescriptor()
    : chemkit::MolecularDescriptor("bpol")
{
    setDimensionality(1);
}

// Returns the sum of the absolute difference between atomic polarizabilities
// of the atoms in each bond in the molecule.
chemkit::Variant BpolDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real value = 0;

    foreach(const chemkit::Bond *bond, molecule->bonds()){
        value += std::abs(atomicPolarizabilities[bond->atom1()->atomicNumber()] -
                          atomicPolarizabilities[bond->atom2()->atomicNumber()]);
    }

    return value;
}
