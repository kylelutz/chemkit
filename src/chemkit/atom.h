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

#ifndef CHEMKIT_ATOM_H
#define CHEMKIT_ATOM_H

#include "chemkit.h"

#include <vector>

#include "point3.h"
#include "element.h"
#include "isotope.h"
#include "vector3.h"
#include "stereochemistry.h"

namespace chemkit {

class Bond;
class Ring;
class Fragment;
class Molecule;
class AtomPrivate;

class CHEMKIT_EXPORT Atom
{
public:
    // properties
    void setElement(const Element &element);
    inline Element element() const;
    void setAtomicNumber(int atomicNumber);
    int atomicNumber() const;
    void setIsotope(const Isotope &isotope);
    Isotope isotope() const;
    void setMassNumber(int massNumber);
    int massNumber() const;
    int expectedValence() const;
    int formalCharge() const;
    void setPartialCharge(Real charge);
    Real partialCharge() const;
    std::string symbol() const;
    std::string name() const;
    Real mass() const;
    Real electronegativity() const;
    Real covalentRadius() const;
    Real vanDerWaalsRadius() const;
    inline bool is(const Element &element) const;
    bool isHeteroatom() const;
    inline Molecule* molecule() const;
    Fragment* fragment() const;
    inline int index() const;

    // structure
    std::vector<Bond *> bonds() const;
    int bondCount() const;
    std::vector<Bond *> bondPathTo(const Atom *atom) const;
    int bondCountTo(const Atom *atom) const;
    int bondCountTo(const Atom *atom, int maxCount) const;
    int valence() const;
    Bond* bondTo(const Atom *atom) const;
    Atom* neighbor(int index) const;
    std::vector<Atom *> neighbors() const;
    int neighborCount() const;
    int neighborCount(const Element &element) const;
    std::vector<Atom *> atomPathTo(const Atom *atom) const;
    int atomCountTo(const Atom *atom) const;
    int atomCountTo(const Atom *atom, int maxCount) const;
    Atom* otherNeighbor(const Atom *neighbor) const;
    bool isBondedTo(const Atom *atom) const;
    bool isBondedTo(const Element &element) const;
    bool isBondedTo(const Element &element, int bondOrder) const;
    bool isConnectedTo(const Atom *atom) const;
    bool isTerminal() const;
    bool isTerminalHydrogen() const;

    // ring perception
    std::vector<Ring *> rings() const;
    int ringCount() const;
    bool isInRing() const;
    bool isInRing(int size) const;
    Ring* smallestRing() const;
    bool isAromatic() const;

    // geometry
    void setPosition(const Point3 &position);
    void setPosition(Real x, Real y, Real z);
    Point3 position() const;
    Real x() const;
    Real y() const;
    Real z() const;
    void moveTo(const Point3 &position);
    void moveTo(Real x, Real y, Real z);
    void moveBy(const Vector3 &vector);
    void moveBy(Real dx, Real dy, Real dz);
    Real distance(const Atom *atom) const;

    // chirality
    void setChirality(Stereochemistry::Type chirality);
    Stereochemistry::Type chirality() const;
    bool isChiral() const;

    enum AtomName{
        Hydrogen = 1,
        Helium = 2,
        Lithium = 3,
        Beryllium = 4,
        Boron = 5,
        Carbon = 6,
        Nitrogen = 7,
        Oxygen = 8,
        Fluorine = 9,
        Neon = 10,
        Sodium = 11,
        Magnesium = 12,
        Aluminum = 13,
        Silicon = 14,
        Phosphorus = 15,
        Sulfur = 16,
        Chlorine = 17,
        Argon = 18,
        Potassium = 19,
        Calcium = 20,
        Scandium = 21,
        Titanium = 22,
        Vanadium = 23,
        Chromium = 24,
        Manganese = 25,
        Iron = 26,
        Cobalt = 27,
        Nickel = 28,
        Copper = 29,
        Zinc = 30,
        Gallium = 31,
        Germanium = 32,
        Arsenic = 33,
        Selenium = 34,
        Bromine = 35,
        Krypton = 36,
        Rubidium = 37,
        Strontium = 38,
        Yttrium = 39,
        Zirconium = 40,
        Niobium = 41,
        Molybdenum = 42,
        Technetium = 43,
        Ruthenium = 44,
        Rhodium = 45,
        Palladium = 46,
        Silver = 47,
        Cadmium = 48,
        Indium = 49,
        Tin = 50,
        Antimony = 51,
        Tellurium = 52,
        Iodine = 53,
        Xenon = 54,
        Cesium = 55,
        Barium = 56,
        Lanthanum = 57,
        Cerium = 58,
        Praseodymium = 59,
        Neodymium = 60,
        Promethium = 61,
        Samarium = 62,
        Europium = 63,
        Gadolinium = 64,
        Terbium = 65,
        Dysprosium = 66,
        Holmium = 67,
        Erbium = 68,
        Thulium = 69,
        Ytterbium = 70,
        Lutetium = 71,
        Hafnium = 72,
        Tantalum = 73,
        Tungsten = 74,
        Rhenium = 75,
        Osmium = 76,
        Iridium = 77,
        Platinum = 78,
        Gold = 79,
        Mercury = 80,
        Thallium = 81,
        Lead = 82,
        Bismuth = 83,
        Polonium = 84,
        Astatine = 85,
        Radon = 86,
        Francium = 87,
        Radium = 88,
        Actinium = 89,
        Thorium = 90,
        Protactinium = 91,
        Uranium = 92,
        Neptunium = 93,
        Plutonium = 94,
        Americium = 95,
        Curium = 96,
        Berkelium = 97,
        Californium = 98,
        Einsteinium = 99,
        Fermium = 100,
        Mendelevium = 101,
        Nobelium = 102,
        Lawrencium = 103,
        Rutherfordium = 104,
        Dubnium = 105,
        Seaborgium = 106,
        Bohrium = 107,
        Hassium = 108,
        Meitnerium = 109
    };

private:
    Atom(Molecule *molecule, const Element &element);
    ~Atom();

    void addBond(Bond *bond);
    void removeBond(Bond *bond);

    CHEMKIT_DISABLE_COPY(Atom)

    friend class Fragment;
    friend class Molecule;

private:
    AtomPrivate* const d;
    Element m_element;
    Molecule *m_molecule;
    int m_index;
};

} // end chemkit namespace

#include "atom-inline.h"

#endif // CHEMKIT_ATOM_H
