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

#include "atommapping.h"

#include "atom.h"
#include "molecule.h"

namespace chemkit {

// === AtomMappingPrivate ================================================== //
class AtomMappingPrivate
{
    public:
        const Molecule *source;
        const Molecule *target;
        QHash<const Atom *, const Atom *> mapping;
};

// === AtomMapping ========================================================= //
/// \class AtomMapping atommapping.h chemkit/atommapping.h
/// \ingroup chemkit
/// \brief The AtomMapping class represents a map between two sets of
///        atoms.

// --- Construction and Destuction ----------------------------------------- //
/// Creates a new, empty atom mapping.
AtomMapping::AtomMapping()
    : d(new AtomMappingPrivate)
{
    d->source = 0;
    d->target = 0;
}

/// Creates a new atom mapping from \p source to \p target.
AtomMapping::AtomMapping(const Molecule *source, const Molecule *target)
    : d(new AtomMappingPrivate)
{
    d->source = source;
    d->target = target;
}

/// Destroys the atom mapping object.
AtomMapping::~AtomMapping()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the source molecule.
const Molecule* AtomMapping::source() const
{
    return d->source;
}

/// Returns the target molecule.
const Molecule* AtomMapping::target() const
{
    return d->target;
}

/// Returns the number of atoms in the mapping.
int AtomMapping::size() const
{
    return d->mapping.size();
}

/// Returns \c true if the mapping is empty (i.e. size() == 0).
bool AtomMapping::isEmpty() const
{
    return size() == 0;
}

// --- Mapping ------------------------------------------------------------- //
/// Adds a new mapping between sourceAtom and targetAtom.
void AtomMapping::add(const Atom *sourceAtom, const Atom *targetAtom)
{
    d->mapping[sourceAtom] = targetAtom;
}

/// Removes the mapping for atom.
void AtomMapping::remove(const Atom *atom)
{
    if(atom->molecule() == d->source){
        d->mapping.remove(atom);
    }
    else if(atom->molecule() == d->target){
        d->mapping.remove(map(atom));
    }
}

/// Returns the atom that atom is mapped to.
const Atom* AtomMapping::map(const Atom *atom) const
{
    if(atom->molecule() == d->source){
        return d->mapping.value(atom, 0);
    }
    else if(atom->molecule() == d->target){
        return d->mapping.key(atom, 0);
    }
    else{
        return 0;
    }
}

/// Removes all atoms in the mapping.
void AtomMapping::clear()
{
    d->mapping.clear();
}

// --- Operators ----------------------------------------------------------- //
AtomMapping& AtomMapping::operator=(const AtomMapping &mapping)
{
    d->source = mapping.d->source;
    d->target = mapping.d->target;
    d->mapping = mapping.d->mapping;

    return *this;
}

} // end chemkit namespace
