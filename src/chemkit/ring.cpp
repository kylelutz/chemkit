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

#include "ring.h"

#include "atom.h"
#include "bond.h"
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
Ring::Ring(QList<Atom *> path)
    : m_atoms(path)
{
    Q_ASSERT(isValid());
}

/// Destroys the ring object.
Ring::~Ring()
{
}

// --- Structure ----------------------------------------------------------- //
/// Returns the number of atoms in the ring with atomicNumber.
int Ring::atomCount(const Element &element) const
{
    int count = 0;

    foreach(const Atom *atom, m_atoms){
        if(atom->is(element)){
            count++;
        }
    }

    return count;
}

/// Returns the bond at \p index in the ring.
Bond* Ring::bond(int index) const
{
    return bonds().value(index);
}

/// Returns the bonds in the ring.
QList<Bond *> Ring::bonds() const
{
    QList<Bond *> bonds;

    for(int i = 0; i < size()-1; i++){
        bonds.append(m_atoms[i]->bondTo(m_atoms[i+1]));
    }

    bonds.append(m_atoms.first()->bondTo(m_atoms.last()));

    return bonds;
}

/// Returns the number of bonds in the ring.
int Ring::bondCount() const
{
    // number of bonds equals the number of atoms
    return m_atoms.size();
}

/// Returns a list of all bonds from atoms inside the ring
/// to atoms outside the ring.
QList<Bond *> Ring::exocyclicBonds() const
{
    QSet<Bond *> bonds;

    foreach(Atom *atom, m_atoms){
        foreach(Bond *bond, atom->bonds()){
            if(!contains(bond)){
                bonds.insert(bond);
            }
        }
    }

    return bonds.toList();
}

/// Returns the number of exocyclic bonds.
int Ring::exocyclicBondCount() const
{
    return exocyclicBonds().size();
}

/// Returns the number of heteroatoms (non-carbon atoms) in the ring.
int Ring::heteroatomCount() const
{
    int count = 0;

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
    QList<Atom *> candidates;

    foreach(Atom *atom, m_atoms){
        if(atom->is(Atom::Carbon))
            continue;

        if(atom->atomicNumber() > highestAtomicNumber){
            candidates.clear();
            candidates.append(atom);
            highestAtomicNumber = atom->atomicNumber();
        }
        else if(atom->atomicNumber() == highestAtomicNumber){
            candidates.append(atom);
        }
    }

    if(candidates.isEmpty()){
        candidates = m_atoms;
    }

    Atom *root = 0;
    int highestNeighborCount = 0;

    foreach(Atom *atom, candidates){
        if(atom->neighborCount() > highestNeighborCount){
            root = atom;
            highestNeighborCount = atom->neighborCount();
        }
    }

    return root;
}

/// Returns the relative postion (distance around the ring) of atom to
/// root. If root is 0, the atom returned from root() is used.
int Ring::position(const Atom *atom, const Atom *root) const
{
    if(!root || !contains(root))
        root = this->root();

    int index = m_atoms.indexOf(const_cast<Atom *>(atom));

    if(index == -1 || atom == root)
        return 0;

    for(int i = 1; i <= (m_atoms.size() / 2); i++){
        int a1, a2;
        a1 = (index + i) % m_atoms.size();
        a2 = (index - i) % m_atoms.size();

        if(a1 < 0) a1 += m_atoms.size();
        if(a2 < 0) a2 += m_atoms.size();

        if(m_atoms[a1] == root || m_atoms[a2] == root)
            return i;
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
    int piCount = 0;

    for(int i = 0; i < size(); i++){
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

    for(int i = 0; i < size(); i++){
        if(!m_atoms[i]->isBondedTo(m_atoms[(i+1) % size()])){
            return false;
        }
    }

    return true;
}

const Atom* Ring::nextAtom(const Atom *atom) const
{
    int index = m_atoms.indexOf(const_cast<Atom *>(atom));

    return m_atoms[(index+1) % size()];
}

const Atom* Ring::previousAtom(const Atom *atom) const
{
    int index = m_atoms.indexOf(const_cast<Atom *>(atom));

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

int Ring::piElectronCount() const
{
    int count = 0;

    for(int i = 0; i < size(); i++){
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
