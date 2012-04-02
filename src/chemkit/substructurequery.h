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

#ifndef CHEMKIT_SUBSTRUCTUREQUERY_H
#define CHEMKIT_SUBSTRUCTUREQUERY_H

#include "chemkit.h"

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "moiety.h"

namespace chemkit {

class Atom;
class Molecule;
class SubstructureQueryPrivate;

class CHEMKIT_EXPORT SubstructureQuery
{
public:
    // enumerations
    enum Flag {
        CompareAtomsOnly = 0x00,
        CompareHydrogens = 0x01,
        CompareAromaticity = 0x02,
        CompareExact = 0x04
    };

    // construction and destruction
    SubstructureQuery();
    SubstructureQuery(const boost::shared_ptr<Molecule> &molecule);
    SubstructureQuery(const std::string &formula, const std::string &format);
    ~SubstructureQuery();

    // properties
    void setMolecule(const boost::shared_ptr<Molecule> &molecule);
    void setMolecule(const std::string &formula, const std::string &format);
    boost::shared_ptr<Molecule> molecule() const;
    void setFlags(int flags);
    int flags() const;

    // queries
    bool matches(const Molecule *molecule) const;
    std::map<Atom *, Atom *> mapping(const Molecule *molecule) const;
    std::map<Atom *, Atom *> maximumMapping(const Molecule *molecule) const;
    std::vector<Molecule *> filter(const std::vector<Molecule *> &molecules) const;
    Moiety find(const Molecule *molecule) const;

private:
    SubstructureQueryPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_SUBSTRUCTUREQUERY_H
