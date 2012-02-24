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

#include "bond.h"

#include <boost/bind.hpp>

#include "atom.h"
#include "ring.h"
#include "foreach.h"
#include "geometry.h"
#include "molecule.h"
#include "moleculeprivate.h"

namespace chemkit {

// === Bond ================================================================ //
/// \class Bond bond.h chemkit/bond.h
/// \ingroup chemkit
/// \brief The Bond class represents a bond between two atoms in a
///        molecule.
///
/// Bond objects are created with the Molecule::addBond() method and
/// destroyed with the Molecule::removeBond() method.

/// \enum Bond::BondType
/// Provides names for the different bond orders:
///     - \c Single = \c 1
///     - \c Double = \c 2
///     - \c Triple = \c 3
///     - \c Quadruple = \c 4

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new bond object for \p molecule at \p index.
Bond::Bond(Molecule *molecule, size_t index)
    : m_molecule(molecule),
      m_index(index)
{
}

/// Destroys the bond object.
Bond::~Bond()
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the atom at \p index in the bond. Index must be either
/// \c 0 or \c 1.
Atom* Bond::atom(size_t index) const
{
    return index == 0 ? atom1() : atom2();
}

/// Returns the first atom in the bond.
Atom* Bond::atom1() const
{
    return m_molecule->d->bondAtoms[m_index].first;
}

/// Returns the second atom in the bond.
Atom* Bond::atom2() const
{
    return m_molecule->d->bondAtoms[m_index].second;
}

/// Returns the other atom in the bond. The given atom must be a
/// part of the bond.
Atom* Bond::otherAtom(const Atom *atom) const
{
    const std::pair<Atom *, Atom *> &pair = m_molecule->d->bondAtoms[m_index];

    return atom == pair.first ? pair.second : pair.first;
}

/// Sets the bond order for the bond to \p order.
void Bond::setOrder(BondOrderType order)
{
    m_molecule->d->bondOrders[m_index] = order;

    molecule()->notifyWatchers(this, MoleculeWatcher::BondOrderChanged);
}

/// Returns the bond order.
Bond::BondOrderType Bond::order() const
{
    return m_molecule->d->bondOrders[m_index];
}

/// Returns \c true if the bond order of the bond is \p order.
bool Bond::is(BondOrderType order) const
{
    return this->order() == order;
}

/// Returns the polarity of the bond. This is calculated as the
/// absolute value of the difference in electronegativity between the
/// two atoms in the bond.
Real Bond::polarity() const
{
    if(atom1()->atomicNumber() == atom2()->atomicNumber()){
        return 0;
    }

    return std::abs(atom1()->electronegativity() - atom2()->electronegativity());
}

/// Returns the dipole moment for the bond.
Vector3 Bond::dipoleMoment() const
{
    Point3 a = atom1()->position();
    Point3 b = atom2()->position();
    Real qa = atom1()->partialCharge();
    Real qb = atom2()->partialCharge();

    return (a - b) * (qa - qb);
}

/// Returns the fragment the bond is a part of.
Fragment* Bond::fragment() const
{
    return atom1()->fragment();
}

// --- Structure ----------------------------------------------------------- //
/// Returns \c true if the bond contains atom.
bool Bond::contains(const Atom *atom) const
{
    return atom1() == atom || atom2() == atom;
}

/// Returns \c true if the bond contains an atom of the given
/// \p element.
bool Bond::contains(const Element &element) const
{
    return atom1()->is(element) || atom2()->is(element);
}

/// Returns \c true if the bond contains both atom \p a and atom
/// \p b.
bool Bond::containsBoth(const Atom *a, const Atom *b) const
{
    return contains(a) && contains(b);
}

/// Returns \c true if the bond contains an atom of the element \p a
/// and an atom of the element \p b.
///
/// For example, to check if this is a carbonyl bond you could use:
/// \code
/// if(bond->containsBoth(Atom::Carbon, Atom::Oxygen) && bond->order() == Bond::Double){
///     // it is a carbonyl
/// }
/// \endcode
bool Bond::containsBoth(const Element &a, const Element &b) const
{
    return (atom1()->is(a) && atom2()->is(b)) || (atom2()->is(a) && atom1()->is(b));
}

/// Returns \c true if either of the two atoms in the bond are
/// terminal.
bool Bond::isTerminal() const
{
    return atom1()->isTerminal() || atom2()->isTerminal();
}

// --- Ring Perception ----------------------------------------------------- //
/// Returns the ring at \p index for the bond.
Ring* Bond::ring(size_t index) const
{
    return rings().advance_begin(index).front();
}

/// Returns a range containing all of the rings that contain the
/// bond.
///
/// \see Molecule::rings()
Bond::RingRange Bond::rings() const
{
    // range of all rings in the molecule
    const Molecule::RingRange moleculeRings = m_molecule->rings();

    // predicate function to check if a ring contains this bond
    boost::function<bool (const Ring *)> ringContainsThisBond =
        boost::bind(static_cast<bool (Ring::*)(const Bond *) const>(&Ring::contains), _1, this);

    return boost::make_iterator_range(
               boost::make_filter_iterator<
                   boost::function<bool (const Ring *)> >(ringContainsThisBond,
                                                          moleculeRings.begin(),
                                                          moleculeRings.end()),
               boost::make_filter_iterator<
                    boost::function<bool (const Ring *)> >(ringContainsThisBond,
                                                           moleculeRings.end(),
                                                           moleculeRings.end()));
}

/// Returns the number of rings that contain the bond.
size_t Bond::ringCount() const
{
    size_t count = 0;

    foreach(const Ring *ring, m_molecule->rings()){
        if(ring->contains(this)){
            count++;
        }
    }

    return count;
}

/// Returns \c true if the bond is a member of at least one ring.
/// (i.e. ringCount() >= 1).
bool Bond::isInRing() const
{
    foreach(const Ring *ring, molecule()->rings()){
        if(ring->contains(this)){
            return true;
        }
    }

    return false;
}

/// Returns \c true if the bond is in a ring of given size.
bool Bond::isInRing(size_t size) const
{
    foreach(const Ring *ring, molecule()->rings()){
        if(ring->size() == size && ring->contains(this)){
            return true;
        }
    }

    return false;
}

/// Returns the smallest ring the bond is a member of or \c 0 if the
/// bond is not in a ring.
Ring* Bond::smallestRing() const
{
    Ring *smallest = 0;

    foreach(Ring *ring, rings()){
        if(!smallest || ring->size() < smallest->size()){
            smallest = ring;
        }
    }

    return smallest;
}

/// Returns \c true if the bond is in an aromatic ring.
///
/// \see Ring::isAromatic()
bool Bond::isAromatic() const
{
    foreach(const Ring *ring, rings()){
        if(ring->isAromatic()){
            return true;
        }
    }

    return false;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the center point of the bond. The center is equal to the
/// midpoint between the two atoms in the bond.
Point3 Bond::center() const
{
    return chemkit::geometry::midpoint(atom1()->position(),
                                       atom2()->position());
}

/// Returns the length of the bond. Length is in Angstroms.
Real Bond::length() const
{
    return atom1()->distance(atom2());
}

// --- Stereochemistry ----------------------------------------------------- //
/// Sets the stereochemistry for the bond.
void Bond::setStereochemistry(Stereochemistry::Type stereochemistry)
{
    m_molecule->stereochemistry()->setStereochemistry(this, stereochemistry);
}

/// Returns the stereochemistry for the bond.
Stereochemistry::Type Bond::stereochemistry() const
{
    if(!m_molecule->m_stereochemistry){
        return Stereochemistry::None;
    }
    else{
        return m_molecule->stereochemistry()->stereochemistry(this);
    }
}

} // end chemkit namespace
