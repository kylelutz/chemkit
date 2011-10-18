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

#include "atom.h"

#include <algorithm>

#include "ring.h"
#include "foreach.h"
#include "vector3.h"
#include "fragment.h"
#include "geometry.h"
#include "molecule.h"
#include "moleculeprivate.h"

namespace chemkit {

// === AtomPrivate ========================================================= //
class AtomPrivate
{
    public:
        Residue *residue;
        Point3 position;
        std::vector<Bond *> bonds;
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
    d->residue = 0;
    d->position = Point3(0, 0, 0);
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
    m_molecule->d->massNumbers[m_index] = massNumber;
    m_molecule->notifyObservers(this, Molecule::AtomMassNumberChanged);
}

/// Returns the mass number of the atom.
int Atom::massNumber() const
{
    return m_molecule->d->massNumbers[m_index];
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
void Atom::setPartialCharge(Real charge)
{
    m_molecule->d->partialCharges[m_index] = charge;
    m_molecule->notifyObservers(this, Molecule::AtomPartialChargeChanged);
}

/// Returns the partial charge of the atom.
Real Atom::partialCharge() const
{
    return m_molecule->d->partialCharges[m_index];
}

/// Returns the elemental symbol for the atom. (e.g. "H" or "Sn").
std::string Atom::symbol() const
{
    return m_element.symbol();
}

/// Returns the elemental name of the atom. (e.g. "Hydrogen" or
/// "Tin").
std::string Atom::name() const
{
    return m_element.name();
}

/// Returns the molar mass of the atom. Mass is in g/mol.
Real Atom::mass() const
{
    return m_element.mass();
}

/// Returns the electronegativity of the atom using the Pauling
/// scale.
Real Atom::electronegativity() const
{
    return m_element.electronegativity();
}

/// Returns the covalent radius of the atom.
Real Atom::covalentRadius() const
{
    return m_element.covalentRadius();
}

/// Returns the Van der Waals radius of the atom.
Real Atom::vanDerWaalsRadius() const
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
    return molecule()->fragment(this);
}

/// Returns the residue the atom is a part of.
Residue* Atom::residue() const
{
    return d->residue;
}

// --- Structure ----------------------------------------------------------- //
/// Returns a list of bonds that this atom is a member of.
std::vector<Bond *> Atom::bonds() const
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
std::vector<Bond *> Atom::bondPathTo(const Atom *atom) const
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
    return neighbors()[index];
}

/// Returns a list of atoms that are directly bonded to the atom.
std::vector<Atom *> Atom::neighbors() const
{
    std::vector<Atom *> neighbors;

    foreach(Bond *bond, d->bonds){
        neighbors.push_back(bond->otherAtom(this));
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
std::vector<Atom *> Atom::atomPathTo(const Atom *atom) const
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
    return fragment()->contains(atom);
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
std::vector<Ring *> Atom::rings() const
{
    std::vector<Ring *> rings;

    foreach(Ring *ring, molecule()->rings()){
        if(ring->contains(this)){
            rings.push_back(ring);
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
void Atom::setPosition(Real x, Real y, Real z)
{
    setPosition(Point3(x, y, z));
}

/// Returns the atom's coordinates.
Point3 Atom::position() const
{
    return d->position;
}

/// Returns the atom's x coordinate. Equivalent to position().x().
Real Atom::x() const
{
    return position().x();
}

/// Returns the atom's y coordinate. Equivalent to position().y().
Real Atom::y() const
{
    return position().y();
}

/// Returns the atom's z coordinate. Equivalent to position().z().
Real Atom::z() const
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
void Atom::moveTo(Real x, Real y, Real z)
{
    setPosition(x, y, z);
}

/// Moves the atom by vector.
void Atom::moveBy(const Vector3 &vector)
{
    moveBy(vector.x(), vector.y(), vector.z());
}

/// Moves the atom by relative amounts dx, dy, dz. Equivalent to
/// setPosition(x()+dx, y()+dy, z()+dz).
void Atom::moveBy(Real dx, Real dy, Real dz)
{
    setPosition(x() + dx, y() + dy, z() + dz);
}

/// Returns the distance between the atom and the other atom.
/// Distance is in Angstroms.
Real Atom::distance(const Atom *atom) const
{
    return chemkit::geometry::distance(position(), atom->position());
}

// --- Chirality ----------------------------------------------------------- //
/// Sets the chirality of the atom.
void Atom::setChirality(Atom::Chirality chirality)
{
    m_molecule->d->chiralities[this] = chirality;
    m_molecule->notifyObservers(this, Molecule::AtomChiralityChanged);
}

/// Returns the chirality of the atom.
Atom::Chirality Atom::chirality() const
{
    const std::map<const Atom*, Atom::Chirality> &chiralities = m_molecule->d->chiralities;

    std::map<const Atom*, Atom::Chirality>::const_iterator location = chiralities.find(this);
    if(location == chiralities.end()){
        return NoChirality;
    }
    else{
        return location->second;
    }
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
    d->bonds.push_back(bond);
}

void Atom::removeBond(Bond *bond)
{
    d->bonds.erase(std::remove(d->bonds.begin(), d->bonds.end(), bond));
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
