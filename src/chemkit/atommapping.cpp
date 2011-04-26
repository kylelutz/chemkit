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
