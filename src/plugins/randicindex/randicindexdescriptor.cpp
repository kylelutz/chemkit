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

#include "randicindexdescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>

RandicIndexDescriptor::RandicIndexDescriptor()
    : chemkit::MolecularDescriptor("randic-index")
{
}

RandicIndexDescriptor::~RandicIndexDescriptor()
{
}

// Returns the randic index for the molecule. See [Randic 1975].
QVariant RandicIndexDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Float value = 0;

    foreach(const chemkit::Bond *bond, molecule->bonds()){
        if(bond->isTerminal() && bond->contains(chemkit::Atom::Hydrogen)){
            continue;
        }

        const chemkit::Atom *a = bond->atom1();
        const chemkit::Atom *b = bond->atom2();

        value += 1.0 / sqrt(chemkit::Float(heavyNeighborCount(a) * heavyNeighborCount(b)));
    }

    return value;
}

// Returns the number of non-hydrogen neighbors for atom.
int RandicIndexDescriptor::heavyNeighborCount(const chemkit::Atom *atom) const
{
    int count = 0;

    foreach(const chemkit::Atom *neighbor, atom->neighbors()){
        if(!neighbor->is(chemkit::Atom::Hydrogen)){
            count++;
        }
    }

    return count;
}
