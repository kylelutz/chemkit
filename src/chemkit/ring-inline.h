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

#ifndef CHEMKIT_RING_INLINE_H
#define CHEMKIT_RING_INLINE_H

#include "ring.h"

#include <algorithm>

#include "atom.h"
#include "bond.h"

namespace chemkit {

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the ring.
inline size_t Ring::size() const
{
    return atomCount();
}

/// Returns the molecule the ring is a part of.
inline Molecule* Ring::molecule() const
{
    return m_atoms[0]->molecule();
}

/// Returns the fragment the ring is a part of.
inline Fragment* Ring::fragment() const
{
    return m_atoms[0]->fragment();
}

// --- Structure ----------------------------------------------------------- //
/// Returns the atom at \p index in the ring.
inline Atom* Ring::atom(size_t index) const
{
    return m_atoms[index];
}

/// Returns the atoms in the ring.
inline std::vector<Atom *> Ring::atoms() const
{
    return m_atoms;
}

/// Returns the number of atoms in the ring.
inline size_t Ring::atomCount() const
{
    return atoms().size();
}

/// Returns an iterator range for the atoms in the ring.
///
/// \internal
inline Ring::AtomRange Ring::atomRange() const
{
    return boost::make_iterator_range(m_atoms.begin(), m_atoms.end());
}

/// Returns \c true if the ring contains atom.
inline bool Ring::contains(const Atom *atom) const
{
    return std::find(m_atoms.begin(), m_atoms.end(), atom) != m_atoms.end();
}

/// Returns \c true if the ring contains bond.
inline bool Ring::contains(const Bond *bond) const
{
    return contains(bond->atom1()) && contains(bond->atom2());
}

/// Returns \c true if the ring contains an atom with atomicNumber.
inline bool Ring::contains(const Element &element) const
{
    for(std::vector<Atom *>::const_iterator iter = m_atoms.begin(); iter != m_atoms.end(); ++iter){
        if((*iter)->is(element)){
            return true;
        }
    }

    return false;
}

} // end chemkit namespace

#endif // CHEMKIT_RING_INLINE_H
