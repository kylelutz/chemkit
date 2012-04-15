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

#ifndef CHEMKIT_MOLECULEPRIVATE_H
#define CHEMKIT_MOLECULEPRIVATE_H

#include "chemkit.h"

#include <map>
#include <string>
#include <vector>

#include "bond.h"
#include "point3.h"
#include "isotope.h"
#include "variantmap.h"

namespace chemkit {

class Atom;
class Bond;
class Ring;
class Fragment;
class CoordinateSet;
class MoleculeWatcher;

class MoleculePrivate
{
public:
    MoleculePrivate();

    std::string name;
    std::vector<Bond *> bonds;
    bool ringsPerceived;
    std::vector<Ring *> rings;
    bool fragmentsPerceived;
    std::vector<Fragment *> fragments;
    std::vector<MoleculeWatcher *> watchers;
    VariantMap data;
    std::map<const Atom *, Isotope> isotopes;
    std::vector<std::string> atomTypes;
    std::vector<Real> partialCharges;
    std::vector<std::pair<Atom*, Atom*> > bondAtoms;
    std::vector<std::vector<Bond *> > atomBonds;
    std::vector<Bond::BondOrderType> bondOrders;
    std::vector<boost::shared_ptr<CoordinateSet> > coordinateSets;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEPRIVATE_H
