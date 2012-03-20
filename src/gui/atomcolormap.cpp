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

#include "atomcolormap.h"

#include <QMap>

#include <chemkit/atom.h>

namespace chemkit {

// === AtomColorMapPrivate ================================================= //
class AtomColorMapPrivate
{
public:
    QMap<int, QColor> colorMap;
    QColor defaultColor;
};

// === AtomColorMap ======================================================== //
/// \class AtomColorMap atomcolormap.h chemkit/atomcolormap.h
/// \ingroup chemkit-gui
/// \brief The AtomColorMap class contains a mapping of elements to
///        colors.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new atom color map.
AtomColorMap::AtomColorMap()
    : d(new AtomColorMapPrivate)
{
}

/// Creates a new atom color map with \p scheme.
AtomColorMap::AtomColorMap(ColorScheme scheme)
    : d(new AtomColorMapPrivate)
{
    setColorScheme(scheme);
}

/// Creates a new atom color map as a copy of \p colorMap.
AtomColorMap::AtomColorMap(const AtomColorMap &colorMap)
    : d(new AtomColorMapPrivate)
{
    d->colorMap = colorMap.d->colorMap;
    d->defaultColor = colorMap.d->defaultColor;
}

/// Destroys the atom color map object.
AtomColorMap::~AtomColorMap()
{
    delete d;
}

// --- Colors -------------------------------------------------------------- //
/// Sets the color for \p element to \p color.
void AtomColorMap::setColor(const Element &element, const QColor &color)
{
    d->colorMap[element.atomicNumber()] = color;
}

/// Returns the color for \p element.
QColor AtomColorMap::color(const Element &element) const
{
    return d->colorMap.value(element.atomicNumber(), defaultColor());
}

/// Returns the color for \p atom.
QColor AtomColorMap::color(const Atom *atom) const
{
    return color(atom->element());
}

/// Sets the default color to \p color.
void AtomColorMap::setDefaultColor(const QColor &color)
{
    d->defaultColor = color;
}

/// Returns the default color for the color map.
QColor AtomColorMap::defaultColor() const
{
    return d->defaultColor;
}

/// Fills the atom color map with the colors according to \p scheme.
void AtomColorMap::setColorScheme(ColorScheme scheme)
{
    // clear existing color map
    d->colorMap.clear();

    if(scheme == DefaultColorScheme){
        d->colorMap[Atom::Hydrogen] = QColor(255, 255, 255);
        d->colorMap[Atom::Helium] = QColor(217, 255, 255);
        d->colorMap[Atom::Lithium] = QColor(204, 128, 255);
        d->colorMap[Atom::Beryllium] = QColor(194, 255, 0);
        d->colorMap[Atom::Boron] = QColor(255, 181, 181);
        d->colorMap[Atom::Carbon] = QColor(80, 80, 80);
        d->colorMap[Atom::Nitrogen] = QColor(48, 80, 248);
        d->colorMap[Atom::Oxygen] = QColor(255, 13, 13);
        d->colorMap[Atom::Fluorine] = QColor(144, 224, 80);
        d->colorMap[Atom::Neon] = QColor(179, 227, 245);
        d->colorMap[Atom::Sodium] = QColor(171, 92, 242);
        d->colorMap[Atom::Magnesium] = QColor(138, 255, 0);
        d->colorMap[Atom::Aluminum] = QColor(191, 166, 166);
        d->colorMap[Atom::Silicon] = QColor(240, 200, 160);
        d->colorMap[Atom::Phosphorus] = QColor(255, 128, 0);
        d->colorMap[Atom::Sulfur] = QColor(255, 255, 48);
        d->colorMap[Atom::Chlorine] = QColor(31, 240, 31);
        d->colorMap[Atom::Argon] = QColor(128, 209, 227);
        d->colorMap[Atom::Potassium] = QColor(143, 64, 212);
        d->colorMap[Atom::Calcium] = QColor(61, 255, 0);
        d->colorMap[Atom::Bromine] = QColor(166, 41, 41);
        d->colorMap[Atom::Iodine] = QColor(148, 0, 148);

        d->defaultColor = QColor(255, 20, 147);
    }
    else if(scheme == RasmolColorScheme){
        d->colorMap[Atom::Carbon] = QColor(200, 200, 200);
        d->colorMap[Atom::Oxygen] = QColor(240, 0, 0);
        d->colorMap[Atom::Hydrogen] = QColor(255, 255, 255);
        d->colorMap[Atom::Nitrogen] = QColor(143, 143, 255);
        d->colorMap[Atom::Sulfur] = QColor(255, 200, 50);
        d->colorMap[Atom::Chlorine] = QColor(0, 255, 0);
        d->colorMap[Atom::Boron] = QColor(0, 255, 0);
        d->colorMap[Atom::Phosphorus] = QColor(255, 165, 0);
        d->colorMap[Atom::Iron] = QColor(255, 165, 0);
        d->colorMap[Atom::Barium] = QColor(255, 165, 0);
        d->colorMap[Atom::Sodium] = QColor(0, 0, 255);
        d->colorMap[Atom::Magnesium] = QColor(34, 139, 34);
        d->colorMap[Atom::Zinc] = QColor(165, 42, 42);
        d->colorMap[Atom::Copper] = QColor(165, 42, 42);
        d->colorMap[Atom::Nickel] = QColor(165, 42, 42);
        d->colorMap[Atom::Bromine] = QColor(165, 42, 42);
        d->colorMap[Atom::Calcium] = QColor(128, 128, 144);
        d->colorMap[Atom::Manganese] = QColor(128, 128, 144);
        d->colorMap[Atom::Aluminum] = QColor(128, 128, 144);
        d->colorMap[Atom::Titanium] = QColor(128, 128, 144);
        d->colorMap[Atom::Chromium] = QColor(128, 128, 144);
        d->colorMap[Atom::Silver] = QColor(128, 128, 144);
        d->colorMap[Atom::Fluorine] = QColor(218, 165, 32);
        d->colorMap[Atom::Silicon] = QColor(218, 165, 32);
        d->colorMap[Atom::Gold] = QColor(218, 165, 32);
        d->colorMap[Atom::Iodine] = QColor(160, 32, 240);
        d->colorMap[Atom::Lithium] = QColor(178, 34, 34);
        d->colorMap[Atom::Helium] = QColor(255, 192, 203);

        d->defaultColor = QColor(255, 20, 147);
    }
    else if(scheme == PymolColorScheme){
        d->colorMap[Atom::Actinium] = QColor(111, 170, 250);
        d->colorMap[Atom::Aluminum] = QColor(191, 165, 165);
        d->colorMap[Atom::Americium] = QColor(84, 92, 242);
        d->colorMap[Atom::Antimony] = QColor(157, 98, 181);
        d->colorMap[Atom::Argon] = QColor(127, 208, 226);
        d->colorMap[Atom::Arsenic] = QColor(189, 127, 226);
        d->colorMap[Atom::Astatine] = QColor(116, 79, 68);
        d->colorMap[Atom::Barium] = QColor(0, 200, 0);
        d->colorMap[Atom::Berkelium] = QColor(138, 79, 226);
        d->colorMap[Atom::Beryllium] = QColor(194, 255, 0);
        d->colorMap[Atom::Bismuth] = QColor(157, 79, 181);
        d->colorMap[Atom::Bohrium] = QColor(224, 0, 55);
        d->colorMap[Atom::Boron] = QColor(255, 181, 181);
        d->colorMap[Atom::Bromine] = QColor(165, 41, 41);
        d->colorMap[Atom::Cadmium] = QColor(255, 216, 143);
        d->colorMap[Atom::Calcium] = QColor(60, 255, 0);
        d->colorMap[Atom::Californium] = QColor(160, 54, 211);
        d->colorMap[Atom::Carbon] = QColor(51, 255, 51);
        d->colorMap[Atom::Cerium] = QColor(255, 255, 199);
        d->colorMap[Atom::Cesium] = QColor(87, 22, 143);
        d->colorMap[Atom::Chlorine] = QColor(30, 240, 30);
        d->colorMap[Atom::Chromium] = QColor(138, 153, 199);
        d->colorMap[Atom::Cobalt] = QColor(240, 143, 159);
        d->colorMap[Atom::Copper] = QColor(199, 127, 51);
        d->colorMap[Atom::Curium] = QColor(119, 92, 226);
        d->colorMap[Atom::Dubnium] = QColor(208, 0, 79);
        d->colorMap[Atom::Dysprosium] = QColor(30, 255, 199);
        d->colorMap[Atom::Einsteinium] = QColor(178, 30, 211);
        d->colorMap[Atom::Erbium] = QColor(0, 229, 116);
        d->colorMap[Atom::Europium] = QColor(97, 255, 199);
        d->colorMap[Atom::Fermium] = QColor(178, 30, 186);
        d->colorMap[Atom::Fluorine] = QColor(178, 255, 255);
        d->colorMap[Atom::Francium] = QColor(65, 0, 102);
        d->colorMap[Atom::Gadolinium] = QColor(68, 255, 199);
        d->colorMap[Atom::Gallium] = QColor(194, 143, 143);
        d->colorMap[Atom::Germanium] = QColor(102, 143, 143);
        d->colorMap[Atom::Gold] = QColor(255, 208, 35);
        d->colorMap[Atom::Hafnium] = QColor(76, 194, 255);
        d->colorMap[Atom::Hassium] = QColor(229, 0, 46);
        d->colorMap[Atom::Helium] = QColor(216, 255, 255);
        d->colorMap[Atom::Holmium] = QColor(0, 255, 156);
        d->colorMap[Atom::Hydrogen] = QColor(229, 229, 229);
        d->colorMap[Atom::Indium] = QColor(165, 116, 114);
        d->colorMap[Atom::Iodine] = QColor(148, 0, 148);
        d->colorMap[Atom::Iridium] = QColor(22, 84, 135);
        d->colorMap[Atom::Iron] = QColor(224, 102, 51);
        d->colorMap[Atom::Krypton] = QColor(92, 183, 208);
        d->colorMap[Atom::Lanthanum] = QColor(111, 211, 255);
        d->colorMap[Atom::Lawrencium] = QColor(199, 0, 102);
        d->colorMap[Atom::Lead] = QColor(87, 89, 97);
        d->colorMap[Atom::Lithium] = QColor(204, 127, 255);
        d->colorMap[Atom::Lutetium] = QColor(0, 170, 36);
        d->colorMap[Atom::Magnesium] = QColor(138, 255, 0);
        d->colorMap[Atom::Manganese] = QColor(156, 122, 199);
        d->colorMap[Atom::Meitnerium] = QColor(234, 0, 38);
        d->colorMap[Atom::Mendelevium] = QColor(178, 12, 165);
        d->colorMap[Atom::Mercury] = QColor(183, 183, 208);
        d->colorMap[Atom::Molybdenum] = QColor(84, 181, 181);
        d->colorMap[Atom::Neodymium] = QColor(199, 255, 199);
        d->colorMap[Atom::Neon] = QColor(178, 226, 245);
        d->colorMap[Atom::Neptunium] = QColor(0, 127, 255);
        d->colorMap[Atom::Nickel] = QColor(79, 208, 79);
        d->colorMap[Atom::Niobium] = QColor(114, 194, 200);
        d->colorMap[Atom::Nitrogen] = QColor(51, 51, 255);
        d->colorMap[Atom::Nobelium] = QColor(189, 12, 135);
        d->colorMap[Atom::Osmium] = QColor(38, 102, 149);
        d->colorMap[Atom::Oxygen] = QColor(255, 76, 76);
        d->colorMap[Atom::Palladium] = QColor(0, 105, 132);
        d->colorMap[Atom::Phosphorus] = QColor(255, 127, 0);
        d->colorMap[Atom::Platinum] = QColor(208, 208, 224);
        d->colorMap[Atom::Plutonium] = QColor(0, 106, 255);
        d->colorMap[Atom::Polonium] = QColor(170, 92, 0);
        d->colorMap[Atom::Potassium] = QColor(143, 63, 211);
        d->colorMap[Atom::Praseodymium] = QColor(216, 255, 199);
        d->colorMap[Atom::Promethium] = QColor(162, 255, 199);
        d->colorMap[Atom::Protactinium] = QColor(0, 160, 255);
        d->colorMap[Atom::Radium] = QColor(0, 124, 0);
        d->colorMap[Atom::Radon] = QColor(65, 130, 149);
        d->colorMap[Atom::Rhenium] = QColor(38, 124, 170);
        d->colorMap[Atom::Rhodium] = QColor(9, 124, 140);
        d->colorMap[Atom::Rubidium] = QColor(111, 46, 175);
        d->colorMap[Atom::Ruthenium] = QColor(36, 143, 143);
        d->colorMap[Atom::Rutherfordium] = QColor(204, 0, 89);
        d->colorMap[Atom::Samarium] = QColor(143, 255, 199);
        d->colorMap[Atom::Scandium] = QColor(229, 229, 229);
        d->colorMap[Atom::Seaborgium] = QColor(216, 0, 68);
        d->colorMap[Atom::Selenium] = QColor(255, 160, 0);
        d->colorMap[Atom::Silicon] = QColor(240, 199, 159);
        d->colorMap[Atom::Silver] = QColor(191, 191, 191);
        d->colorMap[Atom::Sodium] = QColor(170, 92, 242);
        d->colorMap[Atom::Strontium] = QColor(0, 255, 0);
        d->colorMap[Atom::Sulfur] = QColor(229, 197, 63);
        d->colorMap[Atom::Tantalum] = QColor(76, 165, 255);
        d->colorMap[Atom::Technetium] = QColor(58, 157, 157);
        d->colorMap[Atom::Tellurium] = QColor(211, 122, 0);
        d->colorMap[Atom::Terbium] = QColor(47, 255, 199);
        d->colorMap[Atom::Thallium] = QColor(165, 84, 76);
        d->colorMap[Atom::Thorium] = QColor(0, 186, 255);
        d->colorMap[Atom::Thulium] = QColor(0, 211, 81);
        d->colorMap[Atom::Tin] = QColor(102, 127, 127);
        d->colorMap[Atom::Titanium] = QColor(191, 194, 199);
        d->colorMap[Atom::Tungsten] = QColor(33, 148, 213);
        d->colorMap[Atom::Uranium] = QColor(0, 143, 255);
        d->colorMap[Atom::Vanadium] = QColor(165, 165, 170);
        d->colorMap[Atom::Xenon] = QColor(65, 157, 175);
        d->colorMap[Atom::Ytterbium] = QColor(0, 191, 55);
        d->colorMap[Atom::Yttrium] = QColor(148, 255, 255);
        d->colorMap[Atom::Zinc] = QColor(124, 127, 175);
        d->colorMap[Atom::Zirconium] = QColor(148, 224, 224);

        d->defaultColor = QColor(255, 20, 147);
    }
    else if(scheme == JmolColorScheme){
        d->colorMap[Atom::Hydrogen] = QColor(255, 255, 255);
        d->colorMap[Atom::Helium] = QColor(217, 255, 255);
        d->colorMap[Atom::Lithium] = QColor(204, 128, 255);
        d->colorMap[Atom::Beryllium] = QColor(194, 255, 0);
        d->colorMap[Atom::Boron] = QColor(255, 181, 181);
        d->colorMap[Atom::Carbon] = QColor(144, 144, 144);
        d->colorMap[Atom::Nitrogen] = QColor(48, 80, 248);
        d->colorMap[Atom::Oxygen] = QColor(255, 13, 13);
        d->colorMap[Atom::Fluorine] = QColor(144, 224, 80);
        d->colorMap[Atom::Neon] = QColor(179, 227, 245);
        d->colorMap[Atom::Sodium] = QColor(171, 92, 242);
        d->colorMap[Atom::Magnesium] = QColor(138, 255, 0);
        d->colorMap[Atom::Aluminum] = QColor(191, 166, 166);
        d->colorMap[Atom::Silicon] = QColor(240, 200, 160);
        d->colorMap[Atom::Phosphorus] = QColor(255, 128, 0);
        d->colorMap[Atom::Sulfur] = QColor(255, 255, 48);
        d->colorMap[Atom::Chlorine] = QColor(31, 240, 31);
        d->colorMap[Atom::Argon] = QColor(128, 209, 227);
        d->colorMap[Atom::Potassium] = QColor(143, 64, 212);
        d->colorMap[Atom::Calcium] = QColor(61, 255, 0);
        d->colorMap[Atom::Scandium] = QColor(230, 230, 230);
        d->colorMap[Atom::Titanium] = QColor(191, 194, 199);
        d->colorMap[Atom::Vanadium] = QColor(166, 166, 171);
        d->colorMap[Atom::Chromium] = QColor(138, 153, 199);
        d->colorMap[Atom::Manganese] = QColor(156, 122, 199);
        d->colorMap[Atom::Iron] = QColor(224, 102, 51);
        d->colorMap[Atom::Cobalt] = QColor(240, 144, 160);
        d->colorMap[Atom::Nickel] = QColor(80, 208, 80);
        d->colorMap[Atom::Copper] = QColor(200, 128, 51);
        d->colorMap[Atom::Zinc] = QColor(125, 128, 176);
        d->colorMap[Atom::Gallium] = QColor(194, 143, 143);
        d->colorMap[Atom::Germanium] = QColor(102, 143, 143);
        d->colorMap[Atom::Arsenic] = QColor(189, 128, 227);
        d->colorMap[Atom::Selenium] = QColor(255, 161, 0);
        d->colorMap[Atom::Bromine] = QColor(166, 41, 41);
        d->colorMap[Atom::Krypton] = QColor(92, 184, 209);
        d->colorMap[Atom::Rubidium] = QColor(112, 46, 176);
        d->colorMap[Atom::Strontium] = QColor(0, 255, 0);
        d->colorMap[Atom::Yttrium] = QColor(148, 255, 255);
        d->colorMap[Atom::Zirconium] = QColor(148, 224, 224);
        d->colorMap[Atom::Niobium] = QColor(115, 194, 201);
        d->colorMap[Atom::Molybdenum] = QColor(84, 181, 181);
        d->colorMap[Atom::Technetium] = QColor(59, 158, 158);
        d->colorMap[Atom::Ruthenium] = QColor(36, 143, 143);
        d->colorMap[Atom::Rhodium] = QColor(10, 125, 140);
        d->colorMap[Atom::Palladium] = QColor(0, 105, 133);
        d->colorMap[Atom::Silver] = QColor(192, 192, 192);
        d->colorMap[Atom::Cadmium] = QColor(255, 217, 143);
        d->colorMap[Atom::Indium] = QColor(166, 117, 115);
        d->colorMap[Atom::Tin] = QColor(102, 128, 128);
        d->colorMap[Atom::Antimony] = QColor(158, 99, 181);
        d->colorMap[Atom::Tellurium] = QColor(212, 122, 0);
        d->colorMap[Atom::Iodine] = QColor(148, 0, 148);
        d->colorMap[Atom::Xenon] = QColor(66, 158, 176);
        d->colorMap[Atom::Cesium] = QColor(87, 23, 143);
        d->colorMap[Atom::Barium] = QColor(0, 201, 0);
        d->colorMap[Atom::Lanthanum] = QColor(112, 212, 255);
        d->colorMap[Atom::Cerium] = QColor(255, 255, 199);
        d->colorMap[Atom::Praseodymium] = QColor(217, 255, 199);
        d->colorMap[Atom::Neodymium] = QColor(199, 255, 199);
        d->colorMap[Atom::Promethium] = QColor(163, 255, 199);
        d->colorMap[Atom::Samarium] = QColor(143, 255, 199);
        d->colorMap[Atom::Europium] = QColor(97, 255, 199);
        d->colorMap[Atom::Gadolinium] = QColor(69, 255, 199);
        d->colorMap[Atom::Terbium] = QColor(48, 255, 199);
        d->colorMap[Atom::Dysprosium] = QColor(31, 255, 199);
        d->colorMap[Atom::Holmium] = QColor(0, 255, 156);
        d->colorMap[Atom::Erbium] = QColor(0, 230, 117);
        d->colorMap[Atom::Thulium] = QColor(0, 212, 82);
        d->colorMap[Atom::Ytterbium] = QColor(0, 191, 56);
        d->colorMap[Atom::Lutetium] = QColor(0, 171, 36);
        d->colorMap[Atom::Hafnium] = QColor(77, 194, 255);
        d->colorMap[Atom::Tantalum] = QColor(77, 166, 255);
        d->colorMap[Atom::Tungsten] = QColor(33, 148, 214);
        d->colorMap[Atom::Rhenium] = QColor(38, 125, 171);
        d->colorMap[Atom::Osmium] = QColor(38, 102, 150);
        d->colorMap[Atom::Iridium] = QColor(23, 84, 135);
        d->colorMap[Atom::Platinum] = QColor(208, 208, 224);
        d->colorMap[Atom::Gold] = QColor(255, 209, 35);
        d->colorMap[Atom::Mercury] = QColor(184, 184, 208);
        d->colorMap[Atom::Thallium] = QColor(166, 84, 77);
        d->colorMap[Atom::Lead] = QColor(87, 89, 97);
        d->colorMap[Atom::Bismuth] = QColor(158, 79, 181);
        d->colorMap[Atom::Polonium] = QColor(171, 92, 0);
        d->colorMap[Atom::Astatine] = QColor(117, 79, 69);
        d->colorMap[Atom::Radon] = QColor(66, 130, 150);
        d->colorMap[Atom::Francium] = QColor(66, 0, 102);
        d->colorMap[Atom::Radium] = QColor(0, 125, 0);
        d->colorMap[Atom::Actinium] = QColor(112, 171, 250);
        d->colorMap[Atom::Thorium] = QColor(0, 186, 255);
        d->colorMap[Atom::Protactinium] = QColor(0, 161, 255);
        d->colorMap[Atom::Uranium] = QColor(0, 143, 255);
        d->colorMap[Atom::Neptunium] = QColor(0, 128, 255);
        d->colorMap[Atom::Plutonium] = QColor(0, 107, 255);
        d->colorMap[Atom::Americium] = QColor(84, 92, 242);
        d->colorMap[Atom::Curium] = QColor(120, 92, 227);
        d->colorMap[Atom::Berkelium] = QColor(138, 79, 227);
        d->colorMap[Atom::Californium] = QColor(161, 54, 212);
        d->colorMap[Atom::Einsteinium] = QColor(179, 31, 212);
        d->colorMap[Atom::Fermium] = QColor(179, 31, 186);
        d->colorMap[Atom::Mendelevium] = QColor(179, 13, 166);
        d->colorMap[Atom::Nobelium] = QColor(189, 13, 135);
        d->colorMap[Atom::Lawrencium] = QColor(199, 0, 102);
        d->colorMap[Atom::Rutherfordium] = QColor(204, 0, 89);
        d->colorMap[Atom::Dubnium] = QColor(209, 0, 79);
        d->colorMap[Atom::Seaborgium] = QColor(217, 0, 69);
        d->colorMap[Atom::Bohrium] = QColor(224, 0, 56);
        d->colorMap[Atom::Hassium] = QColor(230, 0, 46);
        d->colorMap[Atom::Meitnerium] = QColor(235, 0, 38);

        d->defaultColor = QColor(255, 20, 147);
    }
}

// --- Operators ----------------------------------------------------------- //
AtomColorMap& AtomColorMap::operator=(const AtomColorMap &colorMap)
{
    if(this != &colorMap){
        d->colorMap = colorMap.d->colorMap;
        d->defaultColor = colorMap.d->defaultColor;
    }

    return *this;
}

} // end chemkit namespace
