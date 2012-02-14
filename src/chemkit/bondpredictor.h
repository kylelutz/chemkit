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

#ifndef CHEMKIT_BONDPREDICTOR_H
#define CHEMKIT_BONDPREDICTOR_H

#include "chemkit.h"

#include <vector>

#include <boost/tuple/tuple.hpp>

#include "bond.h"

namespace chemkit {

class Atom;
class Molecule;
class BondPredictorPrivate;

class CHEMKIT_EXPORT BondPredictor
{
public:
    // typedefs
    typedef boost::tuple<Atom *, Atom *, Bond::BondOrderType> PredictedBond;

    // construction and destruction
    BondPredictor(Molecule *molecule);
    ~BondPredictor();

    // properties
    void setTolerance(Real tolerance);
    Real tolerance() const;
    void setMinimumBondLength(Real length);
    Real minimumBondLength() const;
    void setMaximumBondLength(Real length);
    Real maximumBondLength() const;
    Molecule* molecule() const;

    // prediction
    std::vector<PredictedBond> predictedBonds();

    // static methods
    static void predictBonds(Molecule *molecule);

private:
    bool couldBeBonded(Atom *a, Atom *b) const;

private:
    BondPredictorPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_BONDPREDICTOR_H
