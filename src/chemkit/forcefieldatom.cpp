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

#include "forcefieldatom.h"

#include "atom.h"
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
        Float charge;
        Point3 position;
        bool setup;
        ForceField *forceField;
};

// === ForceFieldAtom ====================================================== //
/// \class ForceFieldAtom forcefieldatom.h chemkit/forcefieldatom.h
/// \ingroup chemkit
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
    return forceField()->atoms().indexOf(const_cast<ForceFieldAtom *>(this));
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
void ForceFieldAtom::setCharge(Float charge)
{
    d->charge = charge;
}

/// Returns the charge of the atom.
Float ForceFieldAtom::charge() const
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
Float ForceFieldAtom::energy() const
{
    Float energy = 0;

    Q_FOREACH(const ForceFieldCalculation *calculation, forceField()->calculations()){
        if(calculation->contains(this)){
            energy += calculation->energy();
        }
    }

    return energy;
}

/// Returns the energy gradient for the atom.
Vector3 ForceFieldAtom::gradient() const
{
    return forceField()->gradient()[index()];
}

// --- Structure ----------------------------------------------------------- //
/// Returns \c true if the atom is in a one-four relationship with
/// the other atom.
bool ForceFieldAtom::isOneFour(const ForceFieldAtom *atom) const
{
    const Atom *thisAtom = this->atom();
    const Atom *otherAtom = atom->atom();

    Q_FOREACH(const Atom *neighbor, thisAtom->neighbors()){
        if(neighbor == otherAtom)
            return false;

        Q_FOREACH(const Atom *secondNeighbor, neighbor->neighbors()){
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
    d->position.moveBy(vector);
}

/// Moves the atom's position by (dx, dy, dz).
void ForceFieldAtom::moveBy(Float dx, Float dy, Float dz)
{
    d->position.moveBy(dx, dy, dz);
}

} // end chemkit namespace
