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

#include "element.h"

#include "atom.h"

namespace chemkit {

namespace {

struct ElementData {
    const char *name;
    const char *symbol;
    Real mass; // in grams/mol
    Real electronegativity; // pauling scale
    Real vanDerWaalsRadius; // in angstroms
    Real covalentRadius; // in angstroms
};

const struct ElementData ElementData[] = {
    {"", "", 0, 0, 0, 0},
    {"Hydrogen", "H", 1.00794, 2.1, 1.20, 0.32},
    {"Helium", "He", 4.0026, 0, 1.40, 0.28},
    {"Lithium", "Li", 6.941, 0.98, 1.82, 1.29},
    {"Beryllium", "Be", 9.01218, 1.57, 1.53, 0.96},
    {"Boron", "B", 10.811, 2.04, 1.92, 0.84},
    {"Carbon", "C", 12.011, 2.55, 1.70, 0.77},
    {"Nitrogen", "N", 14.0067, 3.04, 1.55, 0.71},
    {"Oxygen", "O", 15.9994, 3.44, 1.52, 0.66},
    {"Fluorine", "F", 18.9984, 3.98, 1.47, 0.57},
    {"Neon", "Ne", 20.1797, 0, 1.54, 0.58},
    {"Sodium", "Na", 22.98977, 0.93, 2.27, 1.67},
    {"Magnesium", "Mg", 24.305, 1.31, 1.73, 1.42},
    {"Aluminum", "Al", 26.98154, 1.61, 1.84, 1.21},
    {"Silicon", "Si", 28.0855, 1.9, 2.10, 1.11},
    {"Phosphorus", "P", 30.97376, 2.19, 1.80, 1.07},
    {"Sulfur", "S", 32.066, 2.58, 1.80, 1.05},
    {"Chlorine", "Cl", 35.4527, 3.16, 1.75, 1.02},
    {"Argon", "Ar", 39.948, 0, 1.88, 1.06},
    {"Potassium", "K", 39.0983, 0.82, 2.75, 2.03},
    {"Calcium", "Ca", 40.078, 1, 2.31, 1.76},
    {"Scandium", "Sc", 44.9559, 1.36, 2.11, 1.62},
    {"Titanium", "Ti", 47.88, 1.54, 0, 1.32},
    {"Vanadium", "V", 50.9415, 1.63, 0, 1.22},
    {"Chromium", "Cr", 51.996, 1.66, 0, 1.18},
    {"Manganese", "Mn", 54.938, 1.55, 0, 1.17},
    {"Iron", "Fe", 55.847, 1.83, 0, 1.17},
    {"Cobalt", "Co", 58.9332, 1.88, 0, 1.16},
    {"Nickel", "Ni", 58.6934, 1.91, 1.63, 1.15},
    {"Copper", "Cu", 63.546, 1.9, 1.40, 1.28},
    {"Zinc", "Zn", 65.39, 1.65, 1.39, 1.22},
    {"Gallium", "Ga", 69.723, 1.81, 1.87, 1.22},
    {"Germanium", "Ge", 72.61, 2.01, 2.11, 1.22},
    {"Arsenic", "As", 74.9216, 2.18, 1.85, 1.19},
    {"Selenium", "Se", 78.96, 2.55, 1.90, 1.20},
    {"Bromine", "Br", 79.904, 2.96, 1.85, 1.20},
    {"Krypton", "Kr", 83.8, 3, 2.02, 1.16},
    {"Rubidium", "Rb", 85.4678, 0.82, 3.03, 2.20},
    {"Strontium", "Sr", 87.62, 0.95, 2.49, 1.95},
    {"Yttrium", "Y", 88.9059, 1.22, 0, 1.90},
    {"Zirconium", "Zr", 91.224, 1.33, 0, 1.75},
    {"Niobium", "Nb", 92.9064, 1.6, 0, 1.64},
    {"Molybdenum", "Mo", 95.94, 2.16, 0, 1.54},
    {"Technetium", "Tc", 98, 1.9, 0, 1.47},
    {"Ruthenium", "Ru", 101.07, 2.2, 0, 1.46},
    {"Rhodium", "Rh", 102.9055, 2.28, 0, 1.42},
    {"Palladium", "Pd", 106.42, 2.2, 1.63, 1.39},
    {"Silver", "Ag", 107.868, 1.93, 1.72, 1.45},
    {"Cadmium", "Cd", 112.41, 1.69, 1.58, 1.44},
    {"Indium", "In", 114.82, 1.78, 1.93, 1.42},
    {"Tin", "Sn", 118.71, 1.96, 2.17, 1.39},
    {"Antimony", "Sb", 121.757, 2.05, 2.06, 1.39},
    {"Tellurium", "Te", 127.6, 2.1, 2.06, 1.38},
    {"Iodine", "I", 126.9045, 2.66, 1.98, 1.39},
    {"Xenon", "Xe", 131.29, 2.6, 2.16, 1.40},
    {"Cesium", "Cs", 132.9054, 0.79, 3.43, 2.44},
    {"Barium", "Ba", 137.33, 0.89, 2.68, 2.15},
    {"Lanthanum", "La", 138.9055, 1.1, 0, 2.07},
    {"Cerium", "Ce", 140.12, 1.12, 0, 2.04},
    {"Praseodymium", "Pr", 140.9077, 1.13, 0, 2.03},
    {"Neodymium", "Nd", 144.24, 1.14, 0, 2.01},
    {"Promethium", "Pm", 145, 1.13, 0, 1.99},
    {"Samarium", "Sm", 150.36, 1.17, 0, 1.98},
    {"Europium", "Eu", 151.965, 1.2, 0, 1.98},
    {"Gadolinium", "Gd", 157.25, 1.2, 0, 1.96},
    {"Terbium", "Tb", 158.9253, 1.1, 0, 1.94},
    {"Dysprosium", "Dy", 162.5, 1.22, 0, 1.92},
    {"Holmium", "Ho", 164.9303, 1.23, 0, 1.92},
    {"Erbium", "Er", 167.26, 1.24, 0, 1.89},
    {"Thulium", "Tm", 168.9342, 1.25, 0, 1.90},
    {"Ytterbium", "Yb", 173.04, 1.1, 0, 1.87},
    {"Lutetium", "Lu", 174.967, 1.27, 0, 1.76},
    {"Hafnium", "Hf", 178.49, 1.3, 0, 1.75},
    {"Tantalum", "Ta", 180.9479, 1.5, 0, 1.70},
    {"Tungsten", "W", 183.85, 2.36, 0, 1.62},
    {"Rhenium", "Re", 186.207, 1.9, 0, 1.51},
    {"Osmium", "Os", 190.2, 2.2, 0, 1.44},
    {"Iridium", "Ir", 192.22, 2.2, 0, 1.41},
    {"Platinum", "Pt", 195.08, 2.28, 1.75, 1.36},
    {"Gold", "Au", 196.9665, 2.54, 1.66, 1.36},
    {"Mercury", "Hg", 200.59, 2, 1.55, 1.32},
    {"Thallium", "Tl", 204.383, 2.04, 1.96, 1.70},
    {"Lead", "Pb", 207.2, 2.33, 2.02, 1.46},
    {"Bismuth", "Bi", 208.9804, 2.02, 2.07, 1.48},
    {"Polonium", "Po", 209, 2, 1.97, 1.40},
    {"Astatine", "At", 210, 2.2, 2.02, 1.50},
    {"Radon", "Rn", 222, 2.2, 2.20, 1.50},
    {"Francium", "Fr", 223, 0.7, 3.48, 2.60},
    {"Radium", "Ra", 226.0254, 0.9, 2.83, 2.21},
    {"Actinium", "Ac", 227, 1.1, 0, 2.15},
    {"Thorium", "Th", 232.0381, 1.3, 0, 2.06},
    {"Protactinium", "Pa", 231.0359, 1.5, 0, 2.00},
    {"Uranium", "U", 238.029, 1.38, 1.86, 1.96},
    {"Neptunium", "Np", 237.0482, 1.36, 0, 1.90},
    {"Plutonium", "Pu", 244, 1.28, 0, 1.87},
    {"Americium", "Am", 243, 1.3, 0, 1.80},
    {"Curium", "Cm", 247, 1.3, 0, 0},
    {"Berkelium", "Bk", 247, 1.3, 0, 0},
    {"Californium", "Cf", 251, 1.3, 0, 0},
    {"Einsteinium", "Es", 252, 1.3, 0, 0},
    {"Fermium", "Fm", 257, 1.3, 0, 0},
    {"Mendelevium", "Md", 258, 1.3, 0, 0},
    {"Nobelium", "No", 259, 1.3, 0, 0},
    {"Lawrencium", "Lr", 262, 0, 0, 0},
    {"Rutherfordium", "Rf", 261, 0, 0, 0},
    {"Dubnium", "Db", 262, 0, 0, 0},
    {"Seaborgium", "Sg", 263, 0, 0, 0},
    {"Bohrium", "Bh", 262, 0, 0, 0},
    {"Hassium", "Hs", 265, 0, 0, 0},
    {"Meitnerium", "Mt", 266, 0, 0, 0},
};

const int ElementDataSize = sizeof(ElementData) / sizeof(*ElementData);

} // end anonymous namespace

// === Element ============================================================= //
/// \class Element element.h chemkit/element.h
/// \ingroup chemkit
/// \brief The Element class provides information about a chemical
///        element such as its name, symbol, and mass.
///
/// The diagram below shows some of the properties of an element that
/// this class provides access to.
/// \image html element-properties.png "Element Properties" height=2

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new element with the given symbol. If the symbol is not
/// valid the atomic number is set to \c 0.
Element::Element(const char *symbol)
{
    m_atomicNumber = 0;

    for(int i = 1; i < ElementDataSize; i++){
        if(!strcmp(symbol, ElementData[i].symbol)){
            m_atomicNumber = i;
            break;
        }
    }
}

/// Creates a new element with the given symbol. If the symbol is not
/// valid the atomic number is set to \c 0.
Element::Element(const std::string &symbol)
{
    m_atomicNumber = 0;

    for(int i = 1; i < ElementDataSize; i++){
        if(!symbol.compare(ElementData[i].symbol)){
            m_atomicNumber = i;
            break;
        }
    }
}

// --- Properties ---------------------------------------------------------- //
/// Sets the atomic number for the element to \p atomicNumber.
void Element::setAtomicNumber(int atomicNumber)
{
    m_atomicNumber = atomicNumber;
}

/// Returns the element's symbol.
std::string Element::symbol() const
{
    if(!isValid()){
        return std::string();
    }

    return ElementData[m_atomicNumber].symbol;
}

/// Returns the element's name.
std::string Element::name() const
{
    if(!isValid()){
        return std::string();
    }

    return ElementData[m_atomicNumber].name;
}

/// Returns the element's period (row) in the periodic table.
int Element::period() const
{
    if(!isValid())
        return 0;
    else if(m_atomicNumber < 3)
        return 1;
    else if(m_atomicNumber < 11)
        return 2;
    else if(m_atomicNumber < 19)
        return 3;
    else if(m_atomicNumber < 37)
        return 4;
    else if(m_atomicNumber < 55)
        return 5;
    else if(m_atomicNumber < 87)
        return 6;
    else
        return 7;
}

/// Returns the element's mass. Mass is in \c g/mol.
Real Element::mass() const
{
    if(!isValid()){
        return 0;
    }

    return ElementData[m_atomicNumber].mass;
}

/// Returns the element's electronegativity using the Pauling scale.
Real Element::electronegativity() const
{
    if(!isValid()){
        return 0;
    }

    return ElementData[m_atomicNumber].electronegativity;
}

/// Returns the element's covalent radius.
Real Element::covalentRadius() const
{
    if(!isValid()){
        return 0;
    }

    return ElementData[m_atomicNumber].covalentRadius;
}

/// Returns the element's Van der Waals radius.
Real Element::vanDerWaalsRadius() const
{
    if(!isValid()){
        return 0;
    }

    return ElementData[m_atomicNumber].vanDerWaalsRadius;
}

/// Returns the element's expected valence.
int Element::expectedValence() const
{
    switch(m_atomicNumber){
        // hydrogen, alkali metals, halogen group
        case 1:
        case 3:
        case 9:
        case 11:
        case 17:
        case 19:
        case 35:
        case 53:
            return 1;
        // boron group
        case 5:
            return 3;
        // carbon group
        case 6:
        case 14:
            return 4;
        // nitrogen group
        case 7:
        case 15:
        case 33:
            return 3;
        // oxygen group
        case 8:
        case 16:
        case 34:
            return 2;
        default:
            return 0;
    }
}

/// Returns \c true if the element is valid.
bool Element::isValid() const
{
    return isValidAtomicNumber(m_atomicNumber);
}

/// Returns \c true if the element is a metal.
bool Element::isMetal() const
{
    if(!isValid()){
        return false;
    }

    return !isNonmetal();
}

/// Returns \c true if the element is a nonmetal.
bool Element::isNonmetal() const
{
    return m_atomicNumber == Atom::Hydrogen ||
           m_atomicNumber == Atom::Helium ||
           (m_atomicNumber >= Atom::Boron && m_atomicNumber <= Atom::Neon) ||
           (m_atomicNumber >= Atom::Silicon && m_atomicNumber <= Atom::Argon) ||
           (m_atomicNumber >= Atom::Germanium && m_atomicNumber <= Atom::Krypton) ||
           (m_atomicNumber >= Atom::Tellurium && m_atomicNumber <= Atom::Xenon) ||
           (m_atomicNumber >= Atom::Astatine && m_atomicNumber <= Atom::Radon);
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the atomic number for \p symbol. If \p symbol is invalid
/// then \c 0 is returned.
int Element::atomicNumber(const std::string &symbol)
{
    return Element(symbol).atomicNumber();
}

/// Returns the atomic number for \p symbol. If \p symbol is invalid
/// then \c 0 is returned.
int Element::atomicNumber(const char *symbol)
{
    return Element(symbol).atomicNumber();
}

/// Returns the atomic number for \p symbol with \p length. If
/// \p symbol is invalid then \c 0 is returned.
int Element::atomicNumber(const char *symbol, int length)
{
    for(int i = 1; i < ElementDataSize; i++){
        if(!strncmp(symbol, ElementData[i].symbol, length)){
            return i;
        }
    }

    return 0;
}

/// Returns the atomic number for \p symbol. If \p symbol is invalid
/// then \c 0 is returned.
int Element::atomicNumber(char symbol)
{
    switch(symbol){
        case 'H': return 1;
        case 'B': return 5;
        case 'C': return 6;
        case 'N': return 7;
        case 'O': return 8;
        case 'F': return 9;
        case 'P': return 15;
        case 'S': return 16;
        case 'K': return 19;
        case 'V': return 23;
        case 'I': return 53;
        case 'W': return 74;
        case 'U': return 92;
        default: return 0;
    }
}

/// Returns \c true if the atomic number is valid.
bool Element::isValidAtomicNumber(int atomicNumber)
{
    return atomicNumber > 0 && atomicNumber < ElementDataSize;
}

/// Returns \c true if the element symbol is valid.
bool Element::isValidSymbol(const std::string &symbol)
{
    return Element(symbol).isValid();
}

} // end chemkit namespace
