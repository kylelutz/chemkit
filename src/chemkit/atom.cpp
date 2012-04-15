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

#include <boost/bind.hpp>

#include "ring.h"
#include "foreach.h"
#include "vector3.h"
#include "fragment.h"
#include "geometry.h"
#include "molecule.h"
#include "moleculeprivate.h"
#include "cartesiancoordinates.h"

namespace chemkit {

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

// --- Construction and Destruction ---------------------------------------- //
/// Create a new atom object.
Atom::Atom(Molecule *molecule, size_t index)
    : m_molecule(molecule),
      m_index(index)
{
}

/// Destroys the atom object.
Atom::~Atom()
{
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
void Atom::setAtomicNumber(AtomicNumberType atomicNumber)
{
    if(atomicNumber == element().atomicNumber()){
        return;
    }

    if(!Element::isValidAtomicNumber(atomicNumber)){
        return;
    }

    m_molecule->m_elements[m_index].setAtomicNumber(atomicNumber);
    m_molecule->notifyWatchers(this, MoleculeWatcher::AtomElementChanged);
}

/// Returns the atomic number of the atom.
Atom::AtomicNumberType Atom::atomicNumber() const
{
    return element().atomicNumber();
}

/// Set the isotope for the atom to \p isotope.
void Atom::setIsotope(const Isotope &isotope)
{
    if(isotope.element() != element()){
        setElement(isotope.element());
    }

    m_molecule->d->isotopes[this] = isotope;
    m_molecule->notifyWatchers(this, MoleculeWatcher::AtomMassNumberChanged);
}

/// Returns the isotope for the atom.
Isotope Atom::isotope() const
{
    std::map<const Atom *, Isotope>::iterator iter = m_molecule->d->isotopes.find(this);
    if(iter == m_molecule->d->isotopes.end()){
        return Isotope(element(), is(Hydrogen) ? 1 : atomicNumber() * 2);
    }
    else{
        return iter->second;
    }
}

/// Sets the mass number for the atom. This is the number of protons
/// plus the number of neutrons and identifies what isotope the atom
/// is.
void Atom::setMassNumber(MassNumberType massNumber)
{
    setIsotope(Isotope(element(), massNumber));
}

/// Returns the mass number of the atom.
Atom::MassNumberType Atom::massNumber() const
{
    return isotope().massNumber();
}

/// Sets the symbolic type for the atom to \p type.
void Atom::setType(const std::string &type)
{
    std::vector<std::string> &atomTypes = m_molecule->d->atomTypes;

    if(m_index >= atomTypes.size()){
        atomTypes.resize(m_index + 1);
    }

    atomTypes[m_index] = type;
}

/// Returns the symbolic type for the atom or an empty string if no
/// atom type has been set.
std::string Atom::type() const
{
    std::vector<std::string> &atomTypes = m_molecule->d->atomTypes;

    if(m_index >= atomTypes.size()){
        return std::string();
    }

    return atomTypes[m_index];
}

/// Returns the atom's expected valence.
int Atom::expectedValence() const
{
    return element().expectedValence();
}

/// Returns the formal charge of the atom. This is equal to the
/// difference between the atom's valence and its expected valence.
/// (i.e. valence() - expectedValence()).
int Atom::formalCharge() const
{
    if(is(Hydrogen) || element().isMetal())
        return expectedValence() - valence();
    else
        return valence() - expectedValence();
}

/// Sets the partial charge of the atom.
void Atom::setPartialCharge(Real charge)
{
    m_molecule->d->partialCharges[m_index] = charge;
    m_molecule->notifyWatchers(this, MoleculeWatcher::AtomPartialChargeChanged);
}

/// Returns the partial charge of the atom.
Real Atom::partialCharge() const
{
    return m_molecule->d->partialCharges[m_index];
}

/// Returns the elemental symbol for the atom. (e.g. "H" or "Sn").
std::string Atom::symbol() const
{
    return element().symbol();
}

/// Returns the elemental name of the atom. (e.g. "Hydrogen" or
/// "Tin").
std::string Atom::name() const
{
    return element().name();
}

/// Returns the molar mass of the atom. Mass is in g/mol.
Real Atom::mass() const
{
    return element().mass();
}

/// Returns the electronegativity of the atom using the Pauling
/// scale.
Real Atom::electronegativity() const
{
    return element().electronegativity();
}

/// Returns the covalent radius of the atom.
Real Atom::covalentRadius() const
{
    return element().covalentRadius();
}

/// Returns the Van der Waals radius of the atom.
Real Atom::vanDerWaalsRadius() const
{
    return element().vanDerWaalsRadius();
}

/// Returns the fragment the atom is a part of.
Fragment* Atom::fragment() const
{
    return molecule()->fragmentForAtom(this);
}

// --- Structure ----------------------------------------------------------- //
/// Returns the bond at \p index for the atom.
Bond* Atom::bond(size_t index) const
{
    return bonds()[index];
}

/// Returns a range containing all of the bonds that the atom is a
/// member of.
Atom::BondRange Atom::bonds() const
{
    const std::vector<Bond *> &bonds = m_molecule->d->atomBonds[m_index];

    return boost::make_iterator_range(bonds.begin(), bonds.end());
}

/// Returns the number of bonds that this atom is a member of.
/// Equivalent to bonds().size().
size_t Atom::bondCount() const
{
    return bonds().size();
}

/// Returns the number of bonds to the atom.
int Atom::valence() const
{
    int valence = 0;

    foreach(const Bond *bond, bonds()){
        valence += bond->order();
    }

    return valence;
}

/// Returns the bond between the atom and the other atom.
Bond* Atom::bondTo(const Atom *atom) const
{
    foreach(Bond *bond, bonds()){
        if(bond->otherAtom(this) == atom){
            return bond;
        }
    }

    // not bonded to atom
    return 0;
}

/// Returns the bonded neighbor at \p index.
Atom* Atom::neighbor(size_t index) const
{
    return neighbors()[index];
}

/// Returns a range containing all of the atoms that are directly
/// bonded to the atom.
Atom::NeighborRange Atom::neighbors() const
{
    const std::vector<Bond *> &bonds = m_molecule->d->atomBonds[m_index];

    return boost::make_iterator_range(
                boost::make_transform_iterator(
                    bonds.begin(), boost::bind(&Bond::otherAtom, _1, this)),
                boost::make_transform_iterator(
                    bonds.end(), boost::bind(&Bond::otherAtom, _1, this)));
}

/// Returns the number of neighboring (directly bonded) atoms.
size_t Atom::neighborCount() const
{
    return bondCount();
}

/// Returns the number of neighboring atoms of the given \p element.
size_t Atom::neighborCount(const Element &element) const
{
    size_t count = 0;

    foreach(const Bond *bond, bonds()){
        if(bond->otherAtom(this)->is(element)){
            count++;
        }
    }

    return count;
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
/// Returns the ring at \p index for the atom.
Ring* Atom::ring(size_t index) const
{
    return rings().advance_begin(index).front();
}

/// Returns a range containing all of the rings that contain the
/// atom.
///
/// \see Molecule::rings()
Atom::RingRange Atom::rings() const
{
    // range of all rings in the molecule
    const Molecule::RingRange moleculeRings = m_molecule->rings();

    // predicate function to check if a ring contains this atom
    boost::function<bool (const Ring *)> ringContainsThisAtom =
        boost::bind(static_cast<bool (Ring::*)(const Atom *) const>(&Ring::contains), _1, this);

    return boost::make_iterator_range(
               boost::make_filter_iterator<
                   boost::function<bool (const Ring *)> >(ringContainsThisAtom,
                                                          moleculeRings.begin(),
                                                          moleculeRings.end()),
               boost::make_filter_iterator<
                    boost::function<bool (const Ring *)> >(ringContainsThisAtom,
                                                           moleculeRings.end(),
                                                           moleculeRings.end()));
}

/// Returns the number of rings that contain the atom.
size_t Atom::ringCount() const
{
    size_t count = 0;

    foreach(const Ring *ring, m_molecule->rings()){
        if(ring->contains(this)){
            count++;
        }
    }

    return count;
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
bool Atom::isInRing(size_t size) const
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
    m_molecule->coordinates()->setPosition(m_index, position);

    m_molecule->notifyWatchers(this, MoleculeWatcher::AtomPositionChanged);
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
    return m_molecule->coordinates()->position(m_index);
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

/// Returns the distance between the atom and the other atom.
/// Distance is in Angstroms.
Real Atom::distance(const Atom *atom) const
{
    return chemkit::geometry::distance(position(), atom->position());
}

// --- Chirality ----------------------------------------------------------- //
/// Sets the chirality of the atom.
void Atom::setChirality(Stereochemistry::Type chirality)
{
    m_molecule->stereochemistry()->setStereochemistry(this, chirality);
    m_molecule->notifyWatchers(this, MoleculeWatcher::AtomChiralityChanged);
}

/// Returns the chirality of the atom.
Stereochemistry::Type Atom::chirality() const
{
    if(!m_molecule->m_stereochemistry){
        return Stereochemistry::None;
    }
    else{
        return m_molecule->stereochemistry()->stereochemistry(this);
    }
}

/// Returns \c true if the atom is chiral (i.e. chirality() !=
/// Stereochemistry::None).
bool Atom::isChiral() const
{
    return chirality() != Stereochemistry::None;
}

} // end chemkit namespace
