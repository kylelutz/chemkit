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

#include "rotatablebondsdescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>

RotatableBondsDescriptor::RotatableBondsDescriptor()
	: chemkit::MolecularDescriptor("rotatable-bonds")
{
}

RotatableBondsDescriptor::~RotatableBondsDescriptor()
{
}

QVariant RotatableBondsDescriptor::value(const chemkit::Molecule *molecule) const
{
    int count = 0;

    foreach(const chemkit::Bond *bond, molecule->bonds()){
        if(!bond->order() == chemkit::Bond::Single){
            continue;
        }
        else if(heavyNeighborCount(bond->atom1()) < 2 || heavyNeighborCount(bond->atom2()) < 2){
            continue;
        }
        else if(bond->isInRing()){
            continue;
        }

        count++;
    }

    return count;
}

int RotatableBondsDescriptor::heavyNeighborCount(const chemkit::Atom *atom) const
{
    return atom->neighborCount() - atom->neighborCount(chemkit::Atom::Hydrogen);
}
