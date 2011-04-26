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

#include "forcefieldinteractions.h"

#include "foreach.h"
#include "molecule.h"
#include "forcefield.h"

namespace chemkit {

// === ForceFieldInteractionsPrivate ======================================= //
class ForceFieldInteractionsPrivate
{
    public:
        const Molecule *molecule;
        const ForceField *forceField;
};

// === ForceFieldInteractions ============================================== //
/// \class ForceFieldInteractions forcefieldinteractions.h chemkit/forcefieldinteractions.h
/// \ingroup chemkit
/// \internal
/// \brief The ForceFieldInteractions class enumerates atomic
///        interactions in a force field.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new force field interactions object.
ForceFieldInteractions::ForceFieldInteractions(const Molecule *molecule, const ForceField *forceField)
    : d(new ForceFieldInteractionsPrivate)
{
    d->molecule = molecule;
    d->forceField = forceField;
}

/// Destroy the force field interactions object.
ForceFieldInteractions::~ForceFieldInteractions()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the molecule.
const Molecule* ForceFieldInteractions::molecule() const
{
    return d->molecule;
}

/// Returns the force field.
const ForceField* ForceFieldInteractions::forceField() const
{
    return d->forceField;
}

// --- Interactions -------------------------------------------------------- //
/// Returns a list of bonded pairs of atoms.
std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > ForceFieldInteractions::bondedPairs()
{
    std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > bondedPairs;

    Q_FOREACH(const Bond *bond, d->molecule->bonds()){
        const ForceFieldAtom *a = forceField()->atom(bond->atom1());
        const ForceFieldAtom *b = forceField()->atom(bond->atom2());
        bondedPairs.push_back(std::make_pair(a, b));
    }

    return bondedPairs;
}

/// Returns a list of angle groups.
std::vector<std::vector<const ForceFieldAtom *> > ForceFieldInteractions::angleGroups()
{
    std::vector<std::vector<const ForceFieldAtom *> > angleGroups;

    Q_FOREACH(const Atom *atom, d->molecule->atoms()){
        if(!atom->isTerminal()){
            const std::vector<Atom *> &neighbors = atom->neighbors();
            for(unsigned int i = 0; i < neighbors.size(); i++){
                for(unsigned int j = i + 1; j < neighbors.size(); j++){
                    std::vector<const ForceFieldAtom *> angleGroup(3);
                    angleGroup[0] = forceField()->atom(neighbors[i]);
                    angleGroup[1] = forceField()->atom(atom);
                    angleGroup[2] = forceField()->atom(neighbors[j]);

                    angleGroups.push_back(angleGroup);
                }
            }
        }
    }

    return angleGroups;
}

/// Returns a list of torsion groups.
std::vector<std::vector<const ForceFieldAtom *> > ForceFieldInteractions::torsionGroups()
{
    std::vector<std::vector<const ForceFieldAtom *> > torsionGroups;

    std::vector<std::pair<Atom *, Atom *> > pairs;
    Q_FOREACH(const Bond *bond, d->molecule->bonds()){
        if(!bond->atom1()->isTerminal() && !bond->atom2()->isTerminal()){
            pairs.push_back(std::make_pair(bond->atom1(), bond->atom2()));
        }
    }

    // a        d
    //  \      /
    //   b -- c
    std::pair<Atom *, Atom *> pair;
    foreach(pair, pairs){
        const Atom *b = pair.first;
        const Atom *c = pair.second;

        Q_FOREACH(const Atom *a, b->neighbors()){
            if(a == c)
                continue;

            Q_FOREACH(const Atom *d, c->neighbors()){
                if(d == b || d == a)
                    continue;

                std::vector<const ForceFieldAtom *> torsionGroup(4);
                torsionGroup[0] = forceField()->atom(a);
                torsionGroup[1] = forceField()->atom(b);
                torsionGroup[2] = forceField()->atom(c);
                torsionGroup[3] = forceField()->atom(d);

                torsionGroups.push_back(torsionGroup);
            }
        }
    }

    return torsionGroups;
}

/// Returns a list of nonbonded pairs.
std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > ForceFieldInteractions::nonbondedPairs()
{
    std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > nonbondedPairs;

    const std::vector<Atom *> &atoms = d->molecule->atoms();
    for(unsigned int i = 0; i < atoms.size(); i++){
        for(unsigned int j = i + 1; j < atoms.size(); j++){
            if(!atomsWithinTwoBonds(atoms[i], atoms[j])){
                const ForceFieldAtom *a = forceField()->atom(atoms[i]);
                const ForceFieldAtom *b = forceField()->atom(atoms[j]);

                nonbondedPairs.push_back(std::make_pair(a, b));
            }
        }
    }

    return nonbondedPairs;
}

// --- Internal Methods ---------------------------------------------------- //
bool ForceFieldInteractions::atomsWithinTwoBonds(const Atom *a, const Atom *b)
{
    Q_FOREACH(const Atom *neighbor, a->neighbors()){
        if(neighbor == b)
            return true;
        else if(neighbor->isBondedTo(b))
            return true;
    }

    return false;
}

} // end chemkit namespace
