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
QList<QPair<const ForceFieldAtom *, const ForceFieldAtom *> > ForceFieldInteractions::bondedPairs()
{
    QList<QPair<const ForceFieldAtom *, const ForceFieldAtom *> > bondedPairs;

    foreach(const Bond *bond, d->molecule->bonds()){
        const ForceFieldAtom *a = forceField()->atom(bond->atom1());
        const ForceFieldAtom *b = forceField()->atom(bond->atom2());
        bondedPairs.append(qMakePair(a, b));
    }

    return bondedPairs;
}

/// Returns a list of angle groups.
QList<QVector<const ForceFieldAtom *> > ForceFieldInteractions::angleGroups()
{
    QList<QVector<const ForceFieldAtom *> > angleGroups;

    foreach(const Atom *atom, d->molecule->atoms()){
        if(!atom->isTerminal()){
            QList<const Atom *> neighbors = atom->neighbors();
            for(int i = 0; i < neighbors.size(); i++){
                for(int j = i + 1; j < neighbors.size(); j++){
                    QVector<const ForceFieldAtom *> angleGroup(3);
                    angleGroup[0] = forceField()->atom(neighbors[i]);
                    angleGroup[1] = forceField()->atom(atom);
                    angleGroup[2] = forceField()->atom(neighbors[j]);

                    angleGroups.append(angleGroup);
                }
            }
        }
    }

    return angleGroups;
}

/// Returns a list of torsion groups.
QList<QVector<const ForceFieldAtom *> > ForceFieldInteractions::torsionGroups()
{
    QList<QVector<const ForceFieldAtom *> > torsionGroups;

    QList<QPair<const Atom *, const Atom *> > pairs;
    foreach(const Bond *bond, d->molecule->bonds()){
        if(!bond->atom1()->isTerminal() && !bond->atom2()->isTerminal()){
            pairs.append(qMakePair(bond->atom1(), bond->atom2()));
        }
    }

    // a        d
    //  \      /
    //   b -- c
    QPair<const Atom *, const Atom *> pair;
    foreach(pair, pairs){
        const Atom *b = pair.first;
        const Atom *c = pair.second;

        foreach(const Atom *a, b->neighbors()){
            if(a == c)
                continue;

            foreach(const Atom *d, c->neighbors()){
                if(d == b || d == a)
                    continue;

                QVector<const ForceFieldAtom *> torsionGroup(4);
                torsionGroup[0] = forceField()->atom(a);
                torsionGroup[1] = forceField()->atom(b);
                torsionGroup[2] = forceField()->atom(c);
                torsionGroup[3] = forceField()->atom(d);

                torsionGroups.append(torsionGroup);
            }
        }
    }

    return torsionGroups;
}

/// Returns a list of nonbonded pairs.
QList<QPair<const ForceFieldAtom *, const ForceFieldAtom *> > ForceFieldInteractions::nonbondedPairs()
{
    QList<QPair<const ForceFieldAtom *, const ForceFieldAtom *> > nonbondedPairs;

    QList<const Atom *> atoms = d->molecule->atoms();
    for(int i = 0; i < atoms.size(); i++){
        for(int j = i + 1; j < atoms.size(); j++){
            if(!atomsWithinTwoBonds(atoms[i], atoms[j])){
                nonbondedPairs.append(qMakePair(forceField()->atom(atoms[i]),
                                                forceField()->atom(atoms[j])));
            }
        }
    }

    return nonbondedPairs;
}

// --- Internal Methods ---------------------------------------------------- //
bool ForceFieldInteractions::atomsWithinTwoBonds(const Atom *a, const Atom *b)
{
    foreach(const Atom *neighbor, a->neighbors()){
        if(neighbor == b)
            return true;
        else if(neighbor->neighbors().contains(b))
            return true;
    }

    return false;
}

} // end chemkit namespace
