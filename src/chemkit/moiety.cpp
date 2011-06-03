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

#include "moiety.h"
#include "atom.h"

namespace chemkit {

// === MoietyPrivate ======================================================= //
class MoietyPrivate
{
    public:
        std::vector<Atom *> atoms;
};

// === Moiety ============================================================== //
/// \class Moiety moiety.h chemkit/moiety.h
/// \ingroup chemkit
/// \brief The Moiety class represents a group of atoms in a
///        molecule.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty moiety object.
Moiety::Moiety()
    : d(new MoietyPrivate)
{
}

/// Creates a new moiety object containing \p atoms.
Moiety::Moiety(const std::vector<Atom *> &atoms)
    : d(new MoietyPrivate)
{
    d->atoms = atoms;
}

/// Destroys the moiety object.
Moiety::~Moiety()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of atoms in the moiety.
int Moiety::size() const
{
    return atomCount();
}

/// Returns \c true if the moiety contains no atoms.
bool Moiety::isEmpty() const
{
    return size() == 0;
}

/// Returns the molecule that the moiety is a part of. Returns
/// \c 0 if the moiety is empty.
Molecule* Moiety::molecule() const
{
    if(!isEmpty())
        return d->atoms[0]->molecule();
    else
        return 0;

}

// --- Structure ----------------------------------------------------------- //
/// Returns the atom at \p index.
Atom* Moiety::atom(int index) const
{
    return d->atoms[index];
}

/// Returns a list of the atoms in the moiety.
std::vector<Atom *> Moiety::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the moiety.
int Moiety::atomCount() const
{
    return d->atoms.size();
}

// --- Operators ----------------------------------------------------------- //
Moiety& Moiety::operator=(const Moiety &moiety)
{
    d->atoms = moiety.d->atoms;

    return *this;
}

} // end chemkit namespace
