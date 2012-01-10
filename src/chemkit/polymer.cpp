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

#include "polymer.h"

#include "foreach.h"
#include "polymerchain.h"

namespace chemkit {

// === PolymerPrivate ====================================================== //
class PolymerPrivate
{
public:
    std::vector<PolymerChain *> chains;
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
    foreach(PolymerChain *chain, d->chains){
        delete chain;
    }

    delete d;
}

// --- Structure ----------------------------------------------------------- //
/// Adds and returns a new chain to the polymer.
PolymerChain* Polymer::addChain()
{
    PolymerChain *chain = new PolymerChain(this);
    d->chains.push_back(chain);
    return chain;
}

/// Removes the chain from the polymer.
void Polymer::removeChain(PolymerChain *chain)
{
    if(chain->polymer() != this){
        return;
    }

    d->chains.erase(std::remove(d->chains.begin(), d->chains.end(), chain));
    delete chain;
}

/// Returns the chain at \p index.
PolymerChain* Polymer::chain(size_t index) const
{
    return d->chains[index];
}

/// Returns a list of all the chains in the polymer.
std::vector<PolymerChain *> Polymer::chains() const
{
    return d->chains;
}

/// Returns the number of chains in the polymer.
size_t Polymer::chainCount() const
{
    return d->chains.size();
}

} // end chemkit namespace
