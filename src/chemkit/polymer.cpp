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

#include "polymer.h"

#include "polymerchain.h"

namespace chemkit {

// === PolymerPrivate ====================================================== //
class PolymerPrivate
{
    public:
        QList<PolymerChain *> chains;
};

// === Polymer ============================================================= //
/// \class Polymer polymer.h chemkit/polymer.h
/// \ingroup chemkit
/// \brief The Polymer class represents a polymer macromolecule.
///
/// \see PolymerChain

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new polymer.
Polymer::Polymer()
    : Molecule(),
      d(new PolymerPrivate)
{
}

/// Destroys the polymer object.
Polymer::~Polymer()
{
    delete d;
}

// --- Structure ----------------------------------------------------------- //
/// Adds and returns a new chain to the polymer.
PolymerChain* Polymer::addChain()
{
    PolymerChain *chain = new PolymerChain(this);
    d->chains.append(chain);
    return chain;
}

/// Removes the chain from the polymer.
void Polymer::removeChain(PolymerChain *chain)
{
    if(chain->polymer() != this){
        return;
    }

    d->chains.removeOne(chain);
    delete chain;
}

/// Returns the chain at \p index.
PolymerChain* Polymer::chain(int index) const
{
    return d->chains.value(index, 0);
}

/// Returns a list of all the chains in the polymer.
QList<PolymerChain *> Polymer::chains() const
{
    return d->chains;
}

/// Returns the number of chains in the polymer.
int Polymer::chainCount() const
{
    return d->chains.size();
}

} // end chemkit namespace
