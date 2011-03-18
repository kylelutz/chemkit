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

#include "atom.h"

#include "ring.h"
#include "vector.h"
#include "molecule.h"

namespace chemkit {

// === AtomPrivate ========================================================= //
class AtomPrivate
{
    public:
        Residue *residue;
        int massNumber;
        Float partialCharge;
        Point3 position;
        QList<Bond *> bonds;
        Atom::Chirality chirality;
};

// === Atom ================================================================ //
/// \class Atom atom.h chemkit/atom.h
/// \ingroup chemkit
/// \brief The Atom class represents an atom in a molecule.
///
/// Atom objects are created with the Molecule::addAtom() method and destroyed
/// with the Molecule::removeAtom() method.

/// \enum Atom::AtomName
/// Provides names for all the elements. Hydrogen = 1, Helium = 2,
/// Lithium = 3, and so on for all the elements.
///
/// This allows for easier to read code:
/// \code atom->is(Atom::Chlorine); \endcode
/// versus:
/// \code atom->is(17); \endcode

/// \enum Atom::Chirality
/// Provides names for each of the chirality types.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new atom object.
Atom::Atom(Molecule *molecule, const Element &element)
    : d(new AtomPrivate),
      m_element(element),
      m_molecule(molecule)
{
    m_fragment = 0;
    d->residue = 0;

    if(is(Hydrogen)){
        d->massNumber = 1;
    }
    else{
        d->massNumber = element.atomicNumber() * 2;
    }

    d->partialCharge = 0.0;
    d->chirality = NoChirality;
}

/// Destroys the atom object.
Atom::~Atom()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the atom's element to \p element.
void Atom::setElement(const Element &element)
{
    setAtomicNumber(element.atomicNumber());
}

/// Sets the atomic number for the atom. This is the number of
/// protons the atom has and identifies what element the atom is
// (e.g. Hydrogen, Carbon, or Oxygen).
void Atom::setAtomicNumber(int atomicNumber)
{
    if(atomicNumber == m_element.atomicNumber()){
        return;
    }

    if(!Element::isValidAtomicNumber(atomicNumber)){
        return;
    }

    m_element.setAtomicNumber(atomicNumber);
    m_molecule->notifyObservers(this, Molecule::AtomAtomicNumberChanged);
}

/// Returns the atomic number of the atom.
int Atom::atomicNumber() const
{
    return m_element.atomicNumber();
}

/// Sets the mass number for the atom. This is the number of protons
/// plus the number of neutrons and identifies what isotope the atom
/// is.
void Atom::setMassNumber(int massNumber)
{
    d->massNumber = massNumber;
    m_molecule->notifyObservers(this, Molecule::AtomMassNumberChanged);
}

/// Returns the mass number of the atom.
int Atom::massNumber() const
{
    return d->massNumber;
}

/// Returns the atom's expected valence.
int Atom::expectedValence() const
{
    return m_element.expectedValence();
}

/// Returns the formal charge of the atom. This is equal to the
/// difference between the atom's valence and its expected valence.
/// (i.e. valence() - expectedValence()).
int Atom::formalCharge() const
{
    if(is(Hydrogen) || m_element.isMetal())
        return expectedValence() - valence();
    else
        return valence() - expectedValence();
}

/// Sets the partial charge of the atom.
void Atom::setPartialCharge(Float charge)
{
    d->partialCharge = charge;
    m_molecule->notifyObservers(this, Molecule::AtomPartialChargeChanged);
}

/// Returns the partial charge of the atom.
Float Atom::partialCharge() const
{
    return d->partialCharge;
}

/// Returns the elemental symbol for the atom. (e.g. "H" or "Sn").
QString Atom::symbol() const
{
    return m_element.symbol();
}

/// Returns the elemental name of the atom. (e.g. "Hydrogen" or
/// "Tin").
QString Atom::name() const
{
    return m_element.name();
}

/// Returns the molar mass of the atom. Mass is in g/mol.
Float Atom::mass() const
{
    return m_element.mass();
}

/// Returns the electronegativity of the atom using the Pauling
/// scale.
Float Atom::electronegativity() const
{
    return m_element.electronegativity();
}

/// Returns the covalent radius of the atom.
Float Atom::covalentRadius() const
{
    return m_element.covalentRadius();
}

/// Returns the Van der Waals radius of the atom.
Float Atom::vanDerWaalsRadius() const
{
    return m_element.vanDerWaalsRadius();
}

/// Returns \c true if this atom is not Carbon or Hydrogen.
bool Atom::isHeteroatom() const
{
    return !is(Hydrogen) && !is(Carbon);
}

/// Returns the fragment the atom is a part of.
Fragment* Atom::fragment() const
{
    if(m_fragment)
        return m_fragment;
    else
        return molecule()->fragment(this);
}

/// Returns the residue the atom is a part of.
Residue* Atom::residue() const
{
    return d->residue;
}

/// Returns the atom's index.
int Atom::index() const
{
    return m_molecule->indexOf(this);
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of bonds that this atom is a member of.
QList<Bond *> Atom::bonds() const
{
    return d->bonds;
}

/// Returns the number of bonds that this atom is a member of.
/// Equivalent to bonds().size().
int Atom::bondCount() const
{
    return d->bonds.size();
}

/// Returns a list of bonds between the atom and the other atom.
QList<Bond *> Atom::bondPathTo(const Atom *atom) const
{
    return m_molecule->bondPathBetween(this, atom);
}

/// Returns the number of bonds between the atom and the other atom.
int Atom::bondCountTo(const Atom *atom) const
{
    return m_molecule->bondCountBetween(this, atom);
}

/// Returns the number of bonds between the atom and the other atom.
/// Returns \c 0 if \p maxCount is exceeded while searching for
/// \p atom.
int Atom::bondCountTo(const Atom *atom, int maxCount) const
{
    return m_molecule->bondCountBetween(this, atom, maxCount);
}

/// Returns the number of bonds to the atom.
int Atom::valence() const
{
    int valence = 0;

    foreach(const Bond *bond, d->bonds){
        valence += bond->order();
    }

    return valence;
}

/// Returns the bond between the atom and the other atom.
Bond* Atom::bondTo(const Atom *atom) const
{
    foreach(Bond *bond, d->bonds){
        if(bond->otherAtom(this) == atom){
            return bond;
        }
    }

    // not bonded to atom
    return 0;
}

/// Returns the bonded neighbor at \p index.
Atom* Atom::neighbor(int index) const
{
    return neighbors().value(index);
}

/// Returns a list of atoms that are directly bonded to the atom.
QList<Atom *> Atom::neighbors() const
{
    QList<Atom *> neighbors;

    foreach(Bond *bond, d->bonds){
        neighbors.append(bond->otherAtom(this));
    }

    return neighbors;
}

/// Returns the number of neighboring (directly bonded) atoms.
int Atom::neighborCount() const
{
    return bondCount();
}

/// Returns the number of neighboring atoms of the given \p element.
int Atom::neighborCount(const Element &element) const
{
    int count = 0;

    foreach(const Bond *bond, d->bonds){
        if(bond->otherAtom(this)->is(element)){
            count++;
        }
    }

    return count;
}

/// Returns the path of atoms between the atom and the other atom.
QList<Atom *> Atom::atomPathTo(const Atom *atom) const
{
    return m_molecule->atomPathBetween(this, atom);
}

/// Returns the number of atoms between the atom and the other atom.
int Atom::atomCountTo(const Atom *atom) const
{
    return m_molecule->atomCountBetween(this, atom);
}

/// Returns the number of atoms between the atom and the other atom.
/// Returns \c 0 if \p maxCount is exceeded while searching for
/// \p atom.
int Atom::atomCountTo(const Atom *atom, int maxCount) const
{
    return m_molecule->atomCountBetween(this, atom, maxCount);
}

/// Returns the other neighboring atom for a divalent atom.
Atom* Atom::otherNeighbor(const Atom *neighbor) const
{
    foreach(Atom *atom, neighbors()){
        if(atom != neighbor){
            return atom;
        }
    }

    return 0;
}

/// Returns \c true if the atom is bonded to the other atom.
bool Atom::isBondedTo(const Atom *atom) const
{
    return bondTo(atom) != 0;
}

/// Returns \c true if the atom is bonded to an atom of the given
/// \p element.
bool Atom::isBondedTo(const Element &element) const
{
    foreach(const Bond *bond, bonds()){
        if(bond->otherAtom(this)->is(element)){
            return true;
        }
    }

    return false;
}

/// Returns \c true if the atom is bonded to an atom of the given
/// \p element via a bond with \p bondOrder.
bool Atom::isBondedTo(const Element &element, int bondOrder) const
{
    foreach(const Bond *bond, bonds()){
        if(bond->otherAtom(this)->is(element) && bond->order() == bondOrder){
            return true;
        }
    }

    return false;
}

/// Returns \c true if there is a set of contigous bonds that connect
/// this atom with atom (i.e. both atoms are contained in the same
/// fragment).
bool Atom::isConnectedTo(const Atom *atom) const
{
    return atom->fragment() == fragment();
}

/// Returns \c true if this atom is bonded to exactly one atom. (i.e.
/// neighborCount() == 1).
bool Atom::isTerminal() const
{
    return neighborCount() == 1;
}

/// Returns \c true if this atom is bonded to only one atom and is a
/// Hydrogen atom. (i.e. neighborCount() == 1 &&
/// atomicNumber() == 1).
bool Atom::isTerminalHydrogen() const
{
    return isTerminal() && is(Hydrogen);
}

// --- Ring Perception ----------------------------------------------------- //
/// Returns a list of rings the atom is a member of.
///
/// \see Molecule::rings()
QList<Ring *> Atom::rings() const
{
    QList<Ring *> rings;

    foreach(Ring *ring, molecule()->rings()){
        if(ring->contains(this)){
            rings.append(ring);
        }
    }

    return rings;
}

/// Returns the number of rings the atom is a member of.
int Atom::ringCount() const
{
    return rings().size();
}

/// Returns \c true if the atom is a member of at least one ring
/// (i.e. ringCount() >= 1).
bool Atom::isInRing() const
{
    foreach(const Ring *ring, m_molecule->rings()){
        if(ring->contains(this)){
            return true;
        }
    }

    return false;
}

/// Returns \c true if the atom is a member of a ring of given size.
bool Atom::isInRing(int size) const
{
    foreach(const Ring *ring, m_molecule->rings()){
        if(ring->size() == size && ring->contains(this)){
            return true;
        }
    }

    return false;
}

/// Returns the smallest ring the atom is a member of or 0 if the
/// atom is not in a ring.
Ring* Atom::smallestRing() const
{
    Ring *smallest = 0;

    foreach(Ring *ring, rings()){
        if(!smallest || ring->size() < smallest->size()){
            smallest = ring;
        }
    }

    return smallest;
}

/// Returns \c true if the atom is in an aromatic ring.
bool Atom::isAromatic() const
{
    foreach(const Ring *ring, rings()){
        if(ring->isAromatic()){
            return true;
        }
    }

    return false;
}

// --- Geometry ------------------------------------------------------------ //
/// Sets the coordinates of the atom.
void Atom::setPosition(const Point3 &position)
{
    d->position = position;
    molecule()->notifyObservers(this, Molecule::AtomPositionChanged);
}

/// Sets the coordinates of the atom to (x, y, z). Equivalent to
/// setPosition(Point(x, y, z)).
void Atom::setPosition(Float x, Float y, Float z)
{
    setPosition(Point3(x, y, z));
}

/// Returns the atom's coordinates.
Point3 Atom::position() const
{
    return d->position;
}

/// Returns the atom's x coordinate. Equivalent to position().x().
Float Atom::x() const
{
    return position().x();
}

/// Returns the atom's y coordinate. Equivalent to position().y().
Float Atom::y() const
{
    return position().y();
}

/// Returns the atom's z coordinate. Equivalent to position().z().
Float Atom::z() const
{
    return position().z();
}

/// Moves the atom to position. Equivalent to setPosition(position).
void Atom::moveTo(const Point3 &position)
{
    setPosition(position);
}

/// Moves the atom to the point (x, y, z). Equivalent to
/// setPosition(x, y, z).
void Atom::moveTo(Float x, Float y, Float z)
{
    setPosition(x, y, z);
}

/// Moves the atom by vector.
void Atom::moveBy(const Vector &vector)
{
    moveBy(vector.x(), vector.y(), vector.z());
}

/// Moves the atom by relative amounts dx, dy, dz. Equivalent to
/// setPosition(x()+dx, y()+dy, z()+dz).
void Atom::moveBy(Float dx, Float dy, Float dz)
{
    setPosition(x() + dx, y() + dy, z() + dz);
}

/// Returns the distance between the atom and the other atom.
/// Distance is in Angstroms.
Float Atom::distance(const Atom *atom) const
{
    return position().distance(atom->position());
}

// --- Chirality ----------------------------------------------------------- //
/// Sets the chirality of the atom.
void Atom::setChirality(Atom::Chirality chirality)
{
    if(d->chirality == chirality)
        return;

    d->chirality = chirality;
    m_molecule->notifyObservers(this, Molecule::AtomChiralityChanged);
}

/// Returns the chirality of the atom.
Atom::Chirality Atom::chirality() const
{
    return d->chirality;
}

/// Returns \c true if the atom is chiral (i.e. chirality() !=
/// NoChirality).
bool Atom::isChiral() const
{
    return chirality() != NoChirality;
}

// --- Internal Methods ---------------------------------------------------- //
void Atom::addBond(Bond *bond)
{
    Q_ASSERT(d->bonds.contains(bond) == false);

    d->bonds.append(bond);
}

void Atom::removeBond(Bond *bond)
{
    Q_ASSERT(d->bonds.contains(bond) == true);

    d->bonds.removeOne(bond);
}

void Atom::setResidue(Residue *residue)
{
    if(d->residue == residue){
        return;
    }

    d->residue = residue;
    molecule()->notifyObservers(this, Molecule::AtomResidueChanged);
}

} // end chemkit namespace
