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

#include "forcefieldatom.h"

#include <algorithm>

#include <chemkit/atom.h>
#include <chemkit/foreach.h>

#include "forcefield.h"
#include "forcefieldatom.h"
#include "forcefieldcalculation.h"

namespace chemkit {

// === ForceFieldAtomPrivate =============================================== //
class ForceFieldAtomPrivate
{
public:
    const Atom *atom;
    std::string type;
    Real charge;
    Point3 position;
    bool setup;
    ForceField *forceField;
};

// === ForceFieldAtom ====================================================== //
/// \class ForceFieldAtom forcefieldatom.h chemkit/forcefieldatom.h
/// \ingroup chemkit-md
/// \brief The ForceFieldAtom class represents an atom in a force
///        field.
///
/// \see ForceField

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new force field atom.
ForceFieldAtom::ForceFieldAtom(ForceField *forceField, const Atom *atom)
    : d(new ForceFieldAtomPrivate)
{
    d->forceField = forceField;
    d->atom = atom;
    d->position = atom->position();
    d->charge = 0;
    d->setup = false;
}

/// Destroys the force field atom object.
ForceFieldAtom::~ForceFieldAtom()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the atom that the force field atom represents.
const Atom* ForceFieldAtom::atom() const
{
    return d->atom;
}

/// Returns the atom's index.
int ForceFieldAtom::index() const
{
    const std::vector<ForceFieldAtom *> &atoms = forceField()->atoms();

    return std::distance(atoms.begin(), std::find(atoms.begin(), atoms.end(), this));
}

/// Sets the symbolic type for the atom.
bool ForceFieldAtom::setType(const std::string &type)
{
    d->type = type;

    return true;
}

/// Returns the symbolic type for the atom.
std::string ForceFieldAtom::type() const
{
    return d->type;
}

/// Sets the charge of the atom.
void ForceFieldAtom::setCharge(Real charge)
{
    d->charge = charge;
}

/// Returns the charge of the atom.
Real ForceFieldAtom::charge() const
{
    return d->charge;
}

/// Returns \c true if the atom is setup.
bool ForceFieldAtom::isSetup() const
{
    return d->setup;
}

/// Returns the force field the atom is a part of.
ForceField* ForceFieldAtom::forceField() const
{
    return d->forceField;
}

// --- Calculations -------------------------------------------------------- //
/// Returns the total energy of all the calculations the atom is a
/// part of.
Real ForceFieldAtom::energy() const
{
    Real energy = 0;

    foreach(const ForceFieldCalculation *calculation, forceField()->calculations()){
        if(calculation->contains(this)){
            energy += calculation->energy();
        }
    }

    return energy;
}

// --- Structure ----------------------------------------------------------- //
/// Returns \c true if the atom is in a one-four relationship with
/// the other atom.
bool ForceFieldAtom::isOneFour(const ForceFieldAtom *atom) const
{
    const Atom *thisAtom = this->atom();
    const Atom *otherAtom = atom->atom();

    foreach(const Atom *neighbor, thisAtom->neighbors()){
        if(neighbor == otherAtom)
            return false;

        foreach(const Atom *secondNeighbor, neighbor->neighbors()){
            if(secondNeighbor == otherAtom)
                return false;

            if(secondNeighbor->isBondedTo(otherAtom))
                return true;
        }
    }

    return false;
}

// --- Geometry ------------------------------------------------------------ //
/// Sets the position of the atom.
void ForceFieldAtom::setPosition(const Point3 &position)
{
    d->position = position;
}

/// Returns the position of the atom.
Point3 ForceFieldAtom::position() const
{
    return d->position;
}

/// Moves the atom's position by \p vector.
void ForceFieldAtom::moveBy(const Vector3 &vector)
{
    d->position += vector;
}

/// Moves the atom's position by (dx, dy, dz).
void ForceFieldAtom::moveBy(Real dx, Real dy, Real dz)
{
    d->position += Vector3(dx, dy, dz);
}

} // end chemkit namespace
