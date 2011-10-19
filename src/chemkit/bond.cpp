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

#include "atom.h"
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
Bond::Bond(Molecule *molecule, int index)
    : m_molecule(molecule),
      m_index(index)
{
    m_atom1 = 0;
    m_atom2 = 0;
}

/// Destroys the bond object.
Bond::~Bond()
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the atom at \p index in the bond. Index must be either
/// \c 0 or \c 1.
Atom* Bond::atom(int index) const
{
    return index == 0 ? m_atom1 : m_atom2;
}

/// Returns a list containing both atoms in the bond.
std::vector<Atom *> Bond::atoms() const
{
    std::vector<Atom *> atoms(2);
    atoms[0] = m_atom1;
    atoms[1] = m_atom2;
    return atoms;
}

/// Returns the other atom in the bond. The given atom must be a
/// part of the bond.
Atom* Bond::otherAtom(const Atom *atom) const
{
    return m_atom1 == atom ? m_atom2 : m_atom1;
}

/// Sets the bond order for the bond to \p order.
void Bond::setOrder(int order)
{
    m_molecule->d->bondOrders[m_index] = order;

    molecule()->notifyObservers(this, Molecule::BondOrderChanged);
}

/// Returns the bond order.
int Bond::order() const
{
    return m_molecule->d->bondOrders[m_index];
}

/// Returns the polarity of the bond. This is calculated as the
/// absolute value of the difference in electronegativity between the
/// two atoms in the bond.
Real Bond::polarity() const
{
    if(m_atom1->atomicNumber() == m_atom2->atomicNumber()){
        return 0;
    }

    return std::abs(m_atom1->electronegativity() - m_atom2->electronegativity());
}

/// Returns the dipole moment for the bond.
Vector3 Bond::dipoleMoment() const
{
    Point3 a = m_atom1->position();
    Point3 b = m_atom2->position();
    Real qa = m_atom1->partialCharge();
    Real qb = m_atom2->partialCharge();

    return (a - b) * (qa - qb);
}

/// Returns the fragment the bond is a part of.
Fragment* Bond::fragment() const
{
    return m_atom1->fragment();
}

/// Returns the residue the bond is a part of. If the bond is not
/// a part of any residue or the atoms of the bond are in different
/// residues then \c 0 is returned.
Residue* Bond::residue() const
{
    if(m_atom1->residue() == m_atom2->residue())
        return m_atom1->residue();
    else
        return 0;
}

// --- Structure ----------------------------------------------------------- //
/// Returns \c true if the bond contains atom.
bool Bond::contains(const Atom *atom) const
{
    return m_atom1 == atom || m_atom2 == atom;
}

/// Returns \c true if the bond contains an atom of the given
/// \p element.
bool Bond::contains(const Element &element) const
{
    return m_atom1->is(element) || m_atom2->is(element);
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
    return (m_atom1->is(a) && m_atom2->is(b)) || (m_atom2->is(a) && m_atom1->is(b));
}

/// Returns \c true if either of the two atoms in the bond are
/// terminal.
bool Bond::isTerminal() const
{
    return (m_atom1->isTerminal() || m_atom2->isTerminal());
}

// --- Ring Perception ----------------------------------------------------- //
/// Returns a list of rings the bond is a member of.
std::vector<Ring *> Bond::rings() const
{
    std::vector<Ring *> rings;

    foreach(Ring *ring, molecule()->rings()){
        if(ring->contains(this)){
            rings.push_back(ring);
        }
    }

    return rings;
}

/// Returns the number of rings the bond is a member of.
int Bond::ringCount() const
{
    return rings().size();
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
bool Bond::isInRing(int size) const
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
    return chemkit::geometry::midpoint(m_atom1->position(),
                                       m_atom2->position());
}

/// Returns the length of the bond. Length is in Angstroms.
Real Bond::length() const
{
    return m_atom1->distance(m_atom2);
}

} // end chemkit namespace
