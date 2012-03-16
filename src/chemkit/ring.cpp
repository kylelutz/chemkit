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

#include "ring.h"

#include <set>
#include <cassert>
#include <algorithm>

#include <boost/bind.hpp>

#include "atom.h"
#include "bond.h"
#include "foreach.h"
#include "molecule.h"

namespace chemkit {

// === Ring ================================================================ //
/// \class Ring ring.h chemkit/ring.h
/// \ingroup chemkit
/// \brief The Ring class represents a ring of bonded atoms in a
///        molecule.
///
/// Ring objects are returned from the various ring perception
/// methods such as Molecule::rings() and Atom::smallestRing().

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new ring that contains the atoms is \p path.
Ring::Ring(std::vector<Atom *> path)
    : m_atoms(path)
{
    assert(isValid());
}

/// Destroys the ring object.
Ring::~Ring()
{
}

// --- Structure ----------------------------------------------------------- //
/// Returns the number of atoms in the ring with \p element.
size_t Ring::atomCount(const Element &element) const
{
    size_t count = 0;

    foreach(const Atom *atom, m_atoms){
        if(atom->is(element)){
            count++;
        }
    }

    return count;
}

/// Returns the bond at \p index in the ring.
Bond* Ring::bond(size_t index) const
{
    return m_atoms[index]->bondTo(m_atoms[(index + 1) % size()]);
}

/// Returns a range containing all of the bonds in the ring.
Ring::BondRange Ring::bonds() const
{
    return boost::make_iterator_range(
                boost::make_transform_iterator(
                    boost::make_counting_iterator(size_t(0)),
                    boost::bind(&Ring::bond, this, _1)),
                boost::make_transform_iterator(
                    boost::make_counting_iterator(size()),
                    boost::bind(&Ring::bond, this, _1)));
}

/// Returns the number of bonds in the ring.
size_t Ring::bondCount() const
{
    // number of bonds equals the number of atoms
    return m_atoms.size();
}

/// Returns a list of all bonds from atoms inside the ring
/// to atoms outside the ring.
std::vector<Bond *> Ring::exocyclicBonds() const
{
    std::set<Bond *> bondSet;

    foreach(Atom *atom, m_atoms){
        foreach(Bond *bond, atom->bonds()){
            if(!contains(bond)){
                bondSet.insert(bond);
            }
        }
    }

    std::vector<Bond *> bonds;
    foreach(Bond *bond, bondSet){
        bonds.push_back(bond);
    }

    return bonds;
}

/// Returns the number of exocyclic bonds.
size_t Ring::exocyclicBondCount() const
{
    return exocyclicBonds().size();
}

/// Returns the number of heteroatoms (non-carbon atoms) in the ring.
size_t Ring::heteroatomCount() const
{
    size_t count = 0;

    foreach(const Atom *atom, m_atoms){
        if(!atom->is(Atom::Carbon)){
            count++;
        }
    }

    return count;
}

/// Returns \c true if the ring contains any heteroatoms (non-carbon atoms).
bool Ring::isHeterocycle() const
{
    return heteroatomCount() > 0;
}

/// Returns the root atom of the ring. This is determined by finding
/// the non-carbon atom with the highest atomic number. In the case of
/// a tie the atom with the largest valence is returned.
Atom* Ring::root() const
{
    int highestAtomicNumber = 0;
    std::vector<Atom *> candidates;

    foreach(Atom *atom, m_atoms){
        if(atom->is(Atom::Carbon))
            continue;

        if(atom->atomicNumber() > highestAtomicNumber){
            candidates.clear();
            candidates.push_back(atom);
            highestAtomicNumber = atom->atomicNumber();
        }
        else if(atom->atomicNumber() == highestAtomicNumber){
            candidates.push_back(atom);
        }
    }

    if(candidates.empty()){
        candidates = m_atoms;
    }

    Atom *root = 0;
    size_t highestNeighborCount = 0;

    foreach(Atom *atom, candidates){
        if(atom->neighborCount() > highestNeighborCount){
            root = atom;
            highestNeighborCount = atom->neighborCount();
        }
    }

    return root;
}

/// Returns the relative position (distance around the ring) of atom to
/// root. If root is 0, the atom returned from root() is used.
size_t Ring::position(const Atom *atom, const Atom *root) const
{
    if(!root || !contains(root)){
        root = this->root();
    }

    int size = m_atoms.size();
    int index = std::distance(m_atoms.begin(), std::find(m_atoms.begin(), m_atoms.end(), atom));

    if(index == size || atom == root){
        return 0;
    }

    for(int i = 1; i <= (size / 2); i++){
        int a1 = (index + i) % size;
        int a2 = (index - i) % size;

        if(a1 < 0) a1 += size;
        if(a2 < 0) a2 += size;

        if(m_atoms[a1] == root || m_atoms[a2] == root){
            return i;
        }
    }

    // should not get here
    return 0;
}

/// Returns \c true if the ring is fused to \p ring.
bool Ring::isFusedTo(const Ring *ring) const
{
    if(ring == this)
        return false;

    foreach(const Bond *bond, ring->bonds()){
        if(contains(bond)){
            return true;
        }
    }

    return false;
}

// --- Aromaticity --------------------------------------------------------- //
/// Returns \c true if the ring is aromatic.
bool Ring::isAromatic() const
{
    // check for planarity of all ring atoms
    if(!isPlanar()){
        return false;
    }

    // count number of pi electrons
    size_t piCount = 0;

    for(size_t i = 0; i < size(); i++){
        const Atom *atom = m_atoms[i];
        const Bond *nextBond = this->nextBond(atom);

        // double bond
        if(nextBond->order() == Bond::Double)
            piCount += 2;

        // sp2 oxygen group p-orbital lone pair
        if((atom->is(Atom::Oxygen) || atom->is(Atom::Sulfur)) &&
           atom->neighborCount() == 2 && atom->valence() == 2){
            piCount += 2;
        }

        // sp2 nitrogen group p-orbital lone pair
        else if((atom->is(Atom::Nitrogen) || atom->is(Atom::Phosphorus) || atom->is(Atom::Arsenic)) &&
                ((atom->neighborCount() == 3 && atom->valence() == 3) ||
                 (atom->neighborCount() == 2 && atom->valence() == 2))){
            piCount += 2;
        }
    }

    // check for aromaticity using huckel's rule
    if((piCount - 2) % 4 == 0){
        return true;
    }

    // add pi electrons from exocyclic double bonds
    foreach(const Bond *bond, exocyclicBonds()){
        if(bond->order() == Bond::Double){
            piCount += 1;
        }
    }

    // round pi electron count up to next even number if it is odd
    if(piCount & 1){
        piCount++;
    }

    // check again for aromaticity using huckel's rule
    if((piCount - 2) % 4 == 0){
        return true;
    }

    return false;
}

// --- Internal Methods ---------------------------------------------------- //
bool Ring::isValid() const
{
    if(size() < 3)
        return false;

    for(size_t i = 0; i < size(); i++){
        if(!m_atoms[i]->isBondedTo(m_atoms[(i+1) % size()])){
            return false;
        }
    }

    return true;
}

const Atom* Ring::nextAtom(const Atom *atom) const
{
    size_t index = std::distance(m_atoms.begin(), std::find(m_atoms.begin(), m_atoms.end(), atom));

    return m_atoms[(index+1) % size()];
}

const Atom* Ring::previousAtom(const Atom *atom) const
{
    size_t index = std::distance(m_atoms.begin(), std::find(m_atoms.begin(), m_atoms.end(), atom));

    if(index == 0)
        index = size()-1;
    else
        index = (index-1) % size();

    return m_atoms[index];
}

const Bond* Ring::nextBond(const Atom *atom) const
{
    return atom->bondTo(nextAtom(atom));
}

const Bond* Ring::previousBond(const Atom *atom) const
{
    return atom->bondTo(previousAtom(atom));
}

bool Ring::isPlanar() const
{
    foreach(const Atom *atom, m_atoms){
        if(atom->is(Atom::Carbon) && atom->neighborCount() != 3){
            return false;
        }
        else if(atom->is(Atom::Oxygen) && atom->neighborCount() != 2){
            return false;
        }
        else if(atom->is(Atom::Sulfur) && atom->neighborCount() != 2){
            return false;
        }
    }

    return true;
}

size_t Ring::piElectronCount() const
{
    size_t count = 0;

    for(size_t i = 0; i < size(); i++){
        const Atom *atom = m_atoms[i];
        const Bond *nextBond = this->nextBond(atom);
        const Bond *previousBond = this->previousBond(atom);

        // double bond
        if(nextBond->order() == Bond::Double)
            count += 2;

        // sp2 oxygen group p-orbital lone pair
        if((atom->is(Atom::Oxygen) || atom->is(Atom::Sulfur)) &&
           atom->neighborCount() == 2 && atom->valence() == 2){
            count += 2;
        }

        // sp2 nitrogen group p-orbital lone pair
        else if((atom->is(Atom::Nitrogen) || atom->is(Atom::Phosphorus) || atom->is(Atom::Arsenic)) &&
                ((atom->neighborCount() == 3 && atom->valence() == 3) ||
                 (atom->neighborCount() == 2 && atom->valence() == 2))){
            count += 2;
        }

        else{
            // exocyclic double bonds
            foreach(const Bond *bond, atom->bonds()){
                if(bond == nextBond || bond == previousBond){
                    // skip ring bonds
                    continue;
                }

                if(bond->order() == Bond::Double){
                    count += 1;
                }
            }
        }
    }

    return count;
}

} // end chemkit namespace
