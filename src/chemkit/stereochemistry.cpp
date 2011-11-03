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

#include "stereochemistry.h"

#include "atom.h"
#include "bond.h"
#include "molecule.h"

namespace chemkit {

// === Stereochemistry ===================================================== //
/// \class Stereochemistry stereochemistry.h chemkit/stereochemistry.h
/// \ingroup chemkit
/// \brief The Stereochemistry class contains stereochemistry information
///        for the Atom's and Bond's in a Molecule.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new stereochemistry object for \p molecule.
Stereochemistry::Stereochemistry(const Molecule *molecule)
    : m_molecule(molecule)
{
}

/// Destroys the stereochemistry object.
Stereochemistry::~Stereochemistry()
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the molecule for the stereochemistry.
const Molecule* Stereochemistry::molecule() const
{
    return m_molecule;
}

// --- Stereochemistry ----------------------------------------------------- //
/// Sets the stereochemistry for \p atom to \p type.
void Stereochemistry::setStereochemistry(const Atom *atom, Type type)
{
    assert(atom->molecule() == m_molecule);

    if(type == None){
        m_atomStereochemistry.erase(atom);
    }
    else{
        m_atomStereochemistry[atom] = type;
    }
}

/// Sets the stereochemistry for \p bond to \p type.
void Stereochemistry::setStereochemistry(const Bond *bond, Type type)
{
    assert(bond->molecule() == m_molecule);

    if(type == None){
        m_bondStereochemistry.erase(bond);
    }
    else{
        m_bondStereochemistry[bond] = type;
    }
}

/// Returns the stereochemistry for \p atom.
Stereochemistry::Type Stereochemistry::stereochemistry(const Atom *atom) const
{
    std::map<const Atom *, Type>::const_iterator iter = m_atomStereochemistry.find(atom);
    if(iter == m_atomStereochemistry.end()){
        return None;
    }
    else{
        return iter->second;
    }
}

/// Returns the stereochemistry for \p bond.
Stereochemistry::Type Stereochemistry::stereochemistry(const Bond *bond) const
{
    std::map<const Bond *, Type>::const_iterator iter = m_bondStereochemistry.find(bond);
    if(iter == m_bondStereochemistry.end()){
        return None;
    }
    else{
        return iter->second;
    }
}

} // end chemkit namespace
