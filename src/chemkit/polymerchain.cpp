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
    foreach(Residue *residue, d->residues){
        delete residue;
    }

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
size_t PolymerChain::size() const
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
void PolymerChain::insertResidue(size_t index, Residue *residue)
{
    d->residues.insert(d->residues.begin() + index, residue);
}

/// Removes the residue from the chain and deletes it.
bool PolymerChain::removeResidue(Residue *residue)
{
    bool found = takeResidue(residue);

    if(found){
        delete residue;
    }

    return found;
}

/// Removes the residue from the chain.
///
/// The ownership of the residue is passed to the caller.
bool PolymerChain::takeResidue(Residue *residue)
{
    std::vector<Residue *>::iterator location = std::find(d->residues.begin(), d->residues.end(), residue);
    if(location == d->residues.end()){
        return false;
    }

    d->residues.erase(location);
    return true;
}

/// Returns the residue at \p index in the chain.
Residue* PolymerChain::residue(size_t index) const
{
    return d->residues[index];
}

/// Returns a list of the residues in the chain.
std::vector<Residue *> PolymerChain::residues() const
{
    return d->residues;
}

/// Returns the number of residues in the chain.
size_t PolymerChain::residueCount() const
{
    return d->residues.size();
}

/// Returns the index of \p residue in the chain.
int PolymerChain::indexOf(const Residue *residue) const
{
    size_t index = std::distance(d->residues.begin(), std::find(d->residues.begin(), d->residues.end(), residue));
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
