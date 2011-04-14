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
