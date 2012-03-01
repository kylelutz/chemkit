/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#ifndef CHEMKIT_TOPOLOGY_H
#define CHEMKIT_TOPOLOGY_H

#include "md.h"

#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/range/iterator_range.hpp>

namespace chemkit {

class TopologyPrivate;

class CHEMKIT_MD_EXPORT Topology
{
public:
    // typedefs
    typedef boost::array<size_t, 2> BondedInteraction;
    typedef boost::array<size_t, 3> AngleInteraction;
    typedef boost::array<size_t, 4> TorsionInteraction;
    typedef boost::array<size_t, 4> ImproperTorsionInteraction;
    typedef boost::array<size_t, 2> NonbondedInteraction;
    typedef boost::iterator_range<std::vector<BondedInteraction>::const_iterator> BondedInteractionRange;
    typedef boost::iterator_range<std::vector<AngleInteraction>::const_iterator> AngleInteractionRange;
    typedef boost::iterator_range<std::vector<TorsionInteraction>::const_iterator> TorsionInteractionRange;
    typedef boost::iterator_range<std::vector<ImproperTorsionInteraction>::const_iterator> ImproperTorsionInteractionRange;
    typedef boost::iterator_range<std::vector<NonbondedInteraction>::const_iterator> NonbondedInteractionRange;

    // construction and destruction
    Topology();
    Topology(size_t size);
    ~Topology();

    // properties
    void resize(size_t size);
    size_t size() const;
    bool isEmpty() const;

    // atom properties
    void setType(size_t index, const std::string &type);
    std::string type(size_t index) const;
    void setMass(size_t index, Real mass);
    Real mass(size_t index);
    void setCharge(size_t index, Real charge);
    Real charge(size_t index);

    // interations
    void addBondedInteraction(size_t i, size_t j);
    BondedInteractionRange bondedInteractions() const;
    size_t bondedInteractionCount() const;
    void addAngleInteraction(size_t i, size_t j, size_t k);
    AngleInteractionRange angleInteractions() const;
    size_t angleInteractionCount() const;
    void addTorsionInteraction(size_t i, size_t j, size_t k, size_t l);
    TorsionInteractionRange torsionInteractions() const;
    size_t torsionInteractionCount() const;
    void addImproperTorsionInteraction(size_t i, size_t j, size_t k, size_t l);
    ImproperTorsionInteractionRange improperTorsionInteractions() const;
    size_t improperTorsionInteractionCount() const;
    void addNonbondedInteraction(size_t i, size_t j);
    NonbondedInteractionRange nonbondedInteractions() const;
    size_t nonbondedInteractionCount() const;
    bool isOneFour(size_t i, size_t j);

private:
    TopologyPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_TOPOLOGY_H
