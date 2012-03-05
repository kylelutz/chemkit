/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
** Copyright (C) 2005 Chris Morley
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

#include "fp2fingerprint.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>

// The FP2 fingerprint implementation is adapted from code provided
// by Chris Morley.

Fp2Fingerprint::Fp2Fingerprint()
    : chemkit::Fingerprint("fp2")
{
}

Fp2Fingerprint::~Fp2Fingerprint()
{
}

// Returns the FP2 fingerprint value for the molecule.
chemkit::Bitset Fp2Fingerprint::value(const chemkit::Molecule *molecule) const
{
    // create bitset
    chemkit::Bitset fingerprint(1021);

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        // skip fragments starting at terminal hydrogens
        if(atom->isTerminalHydrogen()){
            continue;
        }

        // add each atom fragment to the fingerprint
        addFragments(atom, fingerprint);
    }

    return fingerprint;
}

// Add all fragments starting at atom to the fingerprint.
void Fp2Fingerprint::addFragments(const chemkit::Atom *atom,
                                  chemkit::Bitset &fingerprint) const
{
    Fragment fragment;
    chemkit::Bitset visited(atom->molecule()->atomCount());
    extendFragment(fragment, 1, visited, atom, 0, atom, fingerprint);
}

// Extend the fragment to atom.
void Fp2Fingerprint::extendFragment(Fragment fragment,
                                    size_t depth,
                                    chemkit::Bitset visited,
                                    const chemkit::Atom *atom,
                                    const chemkit::Bond *bond,
                                    const chemkit::Atom *firstAtom,
                                    chemkit::Bitset &fingerprint) const
{
    const size_t MaxFragmentSize = 7;

    chemkit::Bond::BondOrderType bondOrder = 0;
    if(bond){
        bondOrder = bond->isAromatic() ? 5 : bond->order();
    }

    fragment.push_back(bondOrder);
    fragment.push_back(atom->atomicNumber());
    visited.set(atom->index());

    foreach(const chemkit::Bond *neighborBond, atom->bonds()){
        if(neighborBond == bond){
            continue; // don't retrace steps
        }

        const chemkit::Atom *neighbor = neighborBond->otherAtom(atom);
        if(neighbor->isTerminalHydrogen()){
            continue; // don't include terminal hydrogens
        }

        // if the neighbor is an atom that we've already visited
        // then this fragment forms a ring
        if(visited.test(neighbor->index())){
            if(neighbor == firstAtom){
                // add bond at front for the ring
                fragment[0] = bondOrder;

                // insert ring
                Fragment ring = fragment;
                Fragment canonicalRing = fragment;
                for(size_t i = 0; i < ring.size() / 2; i++){
                    // rotate atoms in ring
                    std::rotate(ring.begin(), ring.begin() + 2, ring.end());
                    if(ring > canonicalRing){
                        canonicalRing = ring;
                    }

                    // reverse the ring
                    Fragment reversedRing = ring;
                    std::reverse(reversedRing.begin() + 1, reversedRing.end());
                    if(reversedRing > canonicalRing){
                        canonicalRing = reversedRing;
                    }

                    // add the non-ring form of all ring rotations
                    Fragment ringCopy = ring;
                    ringCopy[0] = 0;
                    fingerprint.set(canonicalHash(ringCopy));
                }

                fingerprint.set(canonicalHash(canonicalRing));
            }
        }
        // no ring
        else{
            if(depth < MaxFragmentSize){
                // extend fragment to the next atom
                extendFragment(fragment,
                               depth+1,
                               visited,
                               neighbor,
                               neighborBond,
                               firstAtom,
                               fingerprint);
            }
        }
    }

    // do not save C, N, O single atom fragments
    if(fragment[0] == 0 && (depth > 1 || fragment[1] > 8 || fragment[1] < 6)){
        fingerprint.set(canonicalHash(fragment));
    }
}

// Returns the canonical hash value for the fragment.
size_t Fp2Fingerprint::canonicalHash(const Fragment &fragment)
{
    const size_t MODINT = 108; // 2^32 % 1021

    // check if we need to reverse the fragment
    bool reverse = false;

    for(size_t i = 1; i < fragment.size(); i++){
        if(*(fragment.begin() + i) == *(fragment.end() - i)){
            continue;
        }

        reverse = *(fragment.begin() + i) < *(fragment.end() - i);
        break;
    }

    // calculate hash value
    size_t hash = 0;

    if(reverse){
        for(size_t i = fragment.size() - 1; i != 0; i--){
            hash = (hash * MODINT + (fragment[i] % 1021)) % 1021;
        }
    }
    else{
        for(size_t i = 0; i < fragment.size(); i++){
            hash = (hash * MODINT + (fragment[i] % 1021)) % 1021;
        }
    }

    return hash;
}
