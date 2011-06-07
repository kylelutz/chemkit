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

#include "conformer.h"

#include <map>

#include "molecule.h"

namespace chemkit {

// === ConformerPrivate ==================================================== //
class ConformerPrivate
{
    public:
        const Molecule *molecule;
        std::map<const Atom *, Point3> coordinates;
};

// === Conformer =========================================================== //
/// \class Conformer conformer.h chemkit/conformer.h
/// \ingroup chemkit
/// \brief The Conformer class represents an alternative set of atomic
///        coordinates for a molecule.
///
/// Conformer objects are created using the Molecule::addConformer()
/// method and destroyed with the Molecule::removeConformer() method.

// --- Construction and Destruction ---------------------------------------- //
Conformer::Conformer(const Molecule *molecule)
    : d(new ConformerPrivate)
{
    d->molecule = molecule;

    Q_FOREACH(const Atom *atom, molecule->atoms()){
        setPosition(atom, atom->position());
    }
}

Conformer::~Conformer()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the molecule for the conformer.
const Molecule* Conformer::molecule() const
{
    return d->molecule;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the coordinates for \p atom to \p position.
void Conformer::setPosition(const Atom *atom, const Point3 &position)
{
    d->coordinates[atom] = position;
}

/// Returns the position of the atom in the conformer.
Point3 Conformer::position(const Atom *atom) const
{
    std::map<const Atom *, Point3>::const_iterator location = d->coordinates.find(atom);
    if(location != d->coordinates.end()){
        return location->second;
    }
    else{
        return atom->position();
    }
}

} // end chemkit namespace
