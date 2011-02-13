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

#ifndef CHEMKIT_BOND_H
#define CHEMKIT_BOND_H

#include "chemkit.h"

#include <QtCore>

#include "point.h"
#include "vector.h"

namespace chemkit {

class Atom;
class Ring;
class Residue;
class Fragment;
class Molecule;

class CHEMKIT_EXPORT Bond
{
    public:
        enum BondType{
            Single = 1,
            Double = 2,
            Triple = 3,
            Quadruple = 4
        };

        // properties
        Atom* atom(int index);
        const Atom* atom(int index) const;
        Atom* atom1();
        const Atom* atom1() const;
        Atom* atom2();
        const Atom* atom2() const;
        QList<Atom *> atoms();
        QList<const Atom *> atoms() const;
        Atom* otherAtom(const Atom *atom);
        const Atom* otherAtom(const Atom *atom) const;
        void setOrder(int order);
        int order() const;
        Float polarity() const;
        Vector dipoleMoment() const;
        Molecule* molecule();
        const Molecule* molecule() const;
        Fragment* fragment();
        const Fragment* fragment() const;
        Residue* residue();
        const Residue* residue() const;
        int index() const;

        // structure
        bool contains(const Atom *atom) const;
        bool contains(int atomicNumber) const;
        bool containsBoth(const Atom *a, const Atom *b) const;
        bool containsBoth(int a, int b) const;
        bool isTerminal() const;

        // ring perception
        QList<Ring *> rings();
        QList<const Ring *> rings() const;
        int ringCount() const;
        bool isInRing() const;
        bool isInRing(int size) const;
        Ring* smallestRing();
        const Ring* smallestRing() const;
        bool isAromatic() const;

        // geometry
        Point center() const;
        Float length() const;

    private:
        Bond(Atom *a, Atom *b, int order = Single);
        ~Bond();

        Q_DISABLE_COPY(Bond)

        friend class Molecule;

    private:
        Atom *m_atom1;
        Atom *m_atom2;
        int m_order;
};

} // end chemkit namespace

#include "bond-inline.h"

#endif // CHEMKIT_BOND_H
