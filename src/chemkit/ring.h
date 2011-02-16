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

#ifndef CHEMKIT_RING_H
#define CHEMKIT_RING_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class Atom;
class Bond;
class Element;
class Fragment;
class Molecule;

class CHEMKIT_EXPORT Ring
{
    public:
        // properties
        int size() const;
        Molecule* molecule() const;
        Fragment* fragment() const;

        // structure
        QList<Atom *> atoms() const;
        int atomCount() const;
        int atomCount(const Element &element) const;
        QList<Bond *> bonds() const;
        int bondCount() const;
        QList<Bond *> exocyclicBonds() const;
        int exocyclicBondCount() const;
        bool contains(const Atom *atom) const;
        bool contains(const Bond *bond) const;
        bool contains(const Element &element) const;
        int heteroatomCount() const;
        bool isHeterocycle() const;
        Atom *root() const;
        int position(const Atom *atom, const Atom *root = 0) const;
        bool isFusedTo(const Ring *ring) const;

        // aromaticity
        bool isAromatic() const;

    private:
        Ring(QList<Atom *> path);
        ~Ring();

        // internal methods
        bool isValid() const;
        const Atom *nextAtom(const Atom *atom) const;
        const Atom *previousAtom(const Atom *atom) const;
        const Bond *nextBond(const Atom *atom) const;
        const Bond *previousBond(const Atom *atom) const;
        bool isPlanar() const;
        int piElectronCount() const;

        Q_DISABLE_COPY(Ring)

        friend class Molecule;
        friend class MolecularGraph;

    private:
        QList<Atom *> m_atoms;
};

} // end chemkit namespace

#include "ring-inline.h"

#endif // CHEMKIT_RING_H
