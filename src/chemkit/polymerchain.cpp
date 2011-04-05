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

#include "polymerchain.h"

#include <sstream>
#include <algorithm>

#include "foreach.h"
#include "residue.h"

namespace chemkit {

// === PolymerChainPrivate ================================================= //
class PolymerChainPrivate
{
    public:
        Polymer *polymer;
        std::string name;
        std::vector<Residue *> residues;
};

// === PolymerChain ======================================================== //
/// \class PolymerChain polymerchain.h chemkit/polymerchain.h
/// \ingroup chemkit
/// \brief The PolymerChain class represents a single chain in a
///        polymer.
///
/// Polymer chains are created with the Polymer::addChain() method
/// and destroyed with the Polymer::removeChain() method.

// --- Construction and Destruction ---------------------------------------- //
PolymerChain::PolymerChain(Polymer *polymer)
    : d(new PolymerChainPrivate)
{
    d->polymer = polymer;
}

PolymerChain::~PolymerChain()
{
    qDeleteAll(d->residues);
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name for the polymer chain to \p name.
void PolymerChain::setName(const std::string &name)
{
    d->name = name;
}

/// Returns the name of the polymer chain.
std::string PolymerChain::name() const
{
    return d->name;
}

/// Returns the number of residues in the chain.
int PolymerChain::size() const
{
    return residueCount();
}

/// Returns \c true if the chain contains no residues.
bool PolymerChain::isEmpty() const
{
    return size() == 0;
}

/// Returns the polymer that the chain belongs to.
Polymer* PolymerChain::polymer() const
{
    return d->polymer;
}

// --- Structure ----------------------------------------------------------- //
/// Adds a residue to the chain.
void PolymerChain::addResidue(Residue *residue)
{
    appendResidue(residue);
}

/// Adds a residue at the end of the chain.
void PolymerChain::appendResidue(Residue *residue)
{
    d->residues.push_back(residue);
}

/// Adds a residue at the beginning of the chain.
void PolymerChain::prependResidue(Residue *residue)
{
    d->residues.insert(d->residues.begin(), residue);
}

/// Adds a residue at \p index of the chain.
///
/// The polymer chain takes ownership of the residue.
void PolymerChain::insertResidue(int index, Residue *residue)
{
    d->residues.insert(d->residues.begin() + index, residue);
}

/// Removes the residue from the chain.
///
/// The ownership of the residue is passed to the caller.
bool PolymerChain::removeResidue(Residue *residue)
{
    std::vector<Residue *>::iterator location = std::find(d->residues.begin(), d->residues.end(), residue);
    if(location == d->residues.end()){
        return false;
    }

    d->residues.erase(location);
    return true;
}

/// Removes the residue from the chain and deletes it.
bool PolymerChain::deleteResidue(Residue *residue)
{
    bool found = removeResidue(residue);

    if(found){
        delete residue;
    }

    return found;
}

/// Returns the residue at \p index in the chain.
Residue* PolymerChain::residue(int index) const
{
    return d->residues[index];
}

/// Returns a list of the residues in the chain.
std::vector<Residue *> PolymerChain::residues() const
{
    return d->residues;
}

/// Returns the number of residues in the chain.
int PolymerChain::residueCount() const
{
    return d->residues.size();
}

/// Returns the index of \p residue in the chain.
int PolymerChain::indexOf(const Residue *residue) const
{
    int index = std::distance(d->residues.begin(), std::find(d->residues.begin(), d->residues.end(), residue));
    if(index == d->residues.size()){
        return -1;
    }

    return index;
}

/// Returns the residue sequence as a string of one letter symbols.
///
/// \see Residue::letter()
std::string PolymerChain::sequenceString() const
{
    std::stringstream stream;

    foreach(const Residue *residue, d->residues){
        stream << residue->letter();
    }

    return stream.str();
}

/// Returns the sequence number of \p residue. Sequence numbers start
/// at \c 1 for the first residue in the chain.
int PolymerChain::sequenceNumber(const Residue *residue) const
{
    int index = indexOf(residue);
    if(index == -1){
        return -1;
    }

    return index + 1;
}

} // end chemkit namespace
