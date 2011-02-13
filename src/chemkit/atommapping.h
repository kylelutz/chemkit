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

#ifndef CHEMKIT_ATOMMAPPING_H
#define CHEMKIT_ATOMMAPPING_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class Atom;
class Molecule;
class AtomMappingPrivate;

class CHEMKIT_EXPORT AtomMapping
{
    public:
        // construction and destruction
        AtomMapping();
        AtomMapping(const Molecule *source, const Molecule *target);
        ~AtomMapping();

        // properties
        const Molecule* source() const;
        const Molecule* target() const;
        int size() const;
        bool isEmpty() const;

        // mapping
        void add(const Atom *sourceAtom, const Atom *targetAtom);
        void remove(const Atom *atom);
        const Atom* map(const Atom *atom) const;
        void clear();

        // operators
        AtomMapping& operator=(const AtomMapping &mapping);

    private:
        AtomMappingPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_ATOMMAPPING_H
