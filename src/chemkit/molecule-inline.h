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

#ifndef CHEMKIT_MOLECULE_INLINE_H
#define CHEMKIT_MOLECULE_INLINE_H

#include "molecule.h"

#include <boost/foreach.hpp>

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the molecule.
inline size_t Molecule::size() const
{
    return atomCount();
}

/// Returns \c true if the molecule contains no atoms (i.e.
/// size() == 0).
inline bool Molecule::isEmpty() const
{
    return size() == 0;
}

// --- Structure ----------------------------------------------------------- //
/// Removes each atom from the molecule which for which \p predicate
/// returns \c true.
template<typename Predicate>
inline void Molecule::removeAtomIf(Predicate predicate)
{
    std::vector<Atom *> atomsToRemove;

    BOOST_FOREACH(Atom *atom, m_atoms){
        if(predicate(atom)){
            atomsToRemove.push_back(atom);
        }
    }

    removeAtoms(atomsToRemove);
}

/// Removes each atom in \p range from the molecule.
template<typename Range>
inline void Molecule::removeAtoms(Range range)
{
    removeAtoms(std::vector<Atom *>(range.begin(), range.end()));
}

/// Returns a range containing all of the atoms in the molecule.
inline Molecule::AtomRange Molecule::atoms() const
{
    return boost::make_iterator_range(m_atoms.begin(), m_atoms.end());
}

/// Returns the number of atoms in the molecule.
inline size_t Molecule::atomCount() const
{
    return m_atoms.size();
}

/// Returns the atom at \p index.
inline Atom* Molecule::atom(size_t index) const
{
    return m_atoms[index];
}

/// Removes each bond from the molecule for which \p predicate
/// returns \c true.
template<typename Predicate>
inline void Molecule::removeBondIf(Predicate predicate)
{
    std::vector<Bond *> bondsToRemove;

    BOOST_FOREACH(Bond *bond, bonds()){
        if(predicate(bond)){
            bondsToRemove.push_back(bond);
        }
    }

    removeBonds(bondsToRemove);
}

/// Removes each bond in \p range from the molecule.
template<typename Range>
inline void Molecule::removeBonds(Range range)
{
    removeBonds(std::vector<Bond *>(range.begin(), range.end()));
}

} // end chemkit namespace

#endif // CHEMKIT_MOLECULE_INLINE_H
