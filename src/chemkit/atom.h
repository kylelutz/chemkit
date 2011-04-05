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

#ifndef CHEMKIT_ATOM_H
#define CHEMKIT_ATOM_H

#include "chemkit.h"

#include <vector>

#include <QtCore>

#include "point3.h"
#include "element.h"
#include "vector3.h"

namespace chemkit {

class Bond;
class Ring;
class Residue;
class Fragment;
class Molecule;
class AtomPrivate;

class CHEMKIT_EXPORT Atom
{
    public:
        // enumerations
        enum Chirality {
            R,
            S,
            NoChirality,
            UnspecifiedChirality
        };

        // properties
        void setElement(const Element &element);
        Element element() const;
        void setAtomicNumber(int atomicNumber);
        int atomicNumber() const;
        void setMassNumber(int massNumber);
        int massNumber() const;
        int expectedValence() const;
        int formalCharge() const;
        void setPartialCharge(Float charge);
        Float partialCharge() const;
        std::string symbol() const;
        std::string name() const;
        Float mass() const;
        Float electronegativity() const;
        Float covalentRadius() const;
        Float vanDerWaalsRadius() const;
        bool is(const Element &element) const;
        bool isHeteroatom() const;
        Molecule* molecule() const;
        Fragment* fragment() const;
        Residue* residue() const;
        int index() const;

        // structure
        std::vector<Bond *> bonds() const;
        int bondCount() const;
        QList<Bond *> bondPathTo(const Atom *atom) const;
        int bondCountTo(const Atom *atom) const;
        int bondCountTo(const Atom *atom, int maxCount) const;
        int valence() const;
        Bond* bondTo(const Atom *atom) const;
        Atom* neighbor(int index) const;
        std::vector<Atom *> neighbors() const;
        int neighborCount() const;
        int neighborCount(const Element &element) const;
        QList<Atom *> atomPathTo(const Atom *atom) const;
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
        void setPosition(Float x, Float y, Float z);
        Point3 position() const;
        Float x() const;
        Float y() const;
        Float z() const;
        void moveTo(const Point3 &position);
        void moveTo(Float x, Float y, Float z);
        void moveBy(const Vector3 &vector);
        void moveBy(Float dx, Float dy, Float dz);
        Float distance(const Atom *atom) const;

        // chirality
        void setChirality(Atom::Chirality chirality);
        Atom::Chirality chirality() const;
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
        void setResidue(Residue *residue);

        Q_DISABLE_COPY(Atom)

        friend class Residue;
        friend class Fragment;
        friend class Molecule;

    private:
        AtomPrivate* const d;
        Element m_element;
        Molecule *m_molecule;
        Fragment *m_fragment;
};

} // end chemkit namespace

#include "atom-inline.h"

#endif // CHEMKIT_ATOM_H
