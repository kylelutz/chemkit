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

#include "topology.h"

#include <chemkit/foreach.h>

namespace chemkit {

// === TopologyPrivate ===================================================== //
class TopologyPrivate
{
public:
    size_t size;
    std::vector<std::string> types;
    std::vector<Real> masses;
    std::vector<Real> charges;
    std::vector<Real> radii;
    std::vector<Topology::BondedInteraction> bondedInteractions;
    std::vector<Topology::AngleInteraction> angleInteractions;
    std::vector<Topology::TorsionInteraction> torsionInteractions;
    std::vector<Topology::ImproperTorsionInteraction> improperTorsionInteractions;
    std::vector<Topology::NonbondedInteraction> nonbondedInteractions;
};

// === Topology ============================================================ //
/// \class Topology topology.h chemkit/topology.h
/// \ingroup chemkit-md
/// \brief The Topology class represents a molecular dynamics topology.
///
/// \see Trajectory, TopologyFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty topology.
Topology::Topology()
    : d(new TopologyPrivate)
{
    d->size = 0;
}

/// Creates a new topology of size \p size.
Topology::Topology(size_t size)
    : d(new TopologyPrivate)
{
    resize(size);
}

/// Destroys the topology object.
Topology::~Topology()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the size of the topolgy to \p size.
void Topology::resize(size_t size)
{
    d->size = size;

    d->types.resize(size);
    d->masses.resize(size);
    d->charges.resize(size);
    d->radii.resize(size);
}

/// Returns the size of the topology.
size_t Topology::size() const
{
    return d->size;
}

/// Returns \c true if the topology is empty (i.e. size() == 0).
bool Topology::isEmpty() const
{
    return size() == 0;
}

// --- Atom Properties ----------------------------------------------------- //
/// Sets the type for the atom at \p index to \p type.
void Topology::setType(size_t index, const std::string &type)
{
    assert(index < d->types.size());

    d->types[index] = type;
}

/// Returns the type for the atom at \p index.
std::string Topology::type(size_t index) const
{
    assert(index < d->types.size());

    return d->types[index];
}

/// Sets the mass for the atom at \p index to \p mass.
void Topology::setMass(size_t index, Real mass)
{
    assert(index < d->masses.size());

    d->masses[index] = mass;
}

/// Returns the mass for the atom at \p index.
Real Topology::mass(size_t index)
{
    assert(index < d->masses.size());

    return d->masses[index];
}

/// Sets the charge for the atom at \p index to \p charge.
void Topology::setCharge(size_t index, Real charge)
{
    assert(index < d->charges.size());

    d->charges[index] = charge;
}

/// Returns the charge for the atom at \p index.
Real Topology::charge(size_t index)
{
    assert(index < d->charges.size());

    return d->charges[index];
}

// --- Interactions -------------------------------------------------------- //
void Topology::addBondedInteraction(size_t i, size_t j)
{
    BondedInteraction interaction;
    interaction[0] = i;
    interaction[1] = j;
    d->bondedInteractions.push_back(interaction);
}

Topology::BondedInteractionRange Topology::bondedInteractions() const
{
    return boost::make_iterator_range(d->bondedInteractions.begin(),
                                      d->bondedInteractions.end());
}

size_t Topology::bondedInteractionCount() const
{
    return d->bondedInteractions.size();
}

void Topology::addAngleInteraction(size_t i, size_t j, size_t k)
{
    AngleInteraction interaction;
    interaction[0] = i;
    interaction[1] = j;
    interaction[2] = k;
    d->angleInteractions.push_back(interaction);
}

Topology::AngleInteractionRange Topology::angleInteractions() const
{
    return boost::make_iterator_range(d->angleInteractions.begin(),
                                      d->angleInteractions.end());
}

size_t Topology::angleInteractionCount() const
{
    return d->angleInteractions.size();
}

void Topology::addTorsionInteraction(size_t i, size_t j, size_t k, size_t l)
{
    TorsionInteraction interaction;
    interaction[0] = i;
    interaction[1] = j;
    interaction[2] = k;
    interaction[3] = l;
    d->torsionInteractions.push_back(interaction);
}

Topology::TorsionInteractionRange Topology::torsionInteractions() const
{
    return boost::make_iterator_range(d->torsionInteractions.begin(),
                                      d->torsionInteractions.end());
}

size_t Topology::torsionInteractionCount() const
{
    return d->torsionInteractions.size();
}

void Topology::addImproperTorsionInteraction(size_t i, size_t j, size_t k, size_t l)
{
    ImproperTorsionInteraction interaction;
    interaction[0] = i;
    interaction[1] = j;
    interaction[2] = k;
    interaction[3] = l;
    d->improperTorsionInteractions.push_back(interaction);
}

Topology::ImproperTorsionInteractionRange Topology::improperTorsionInteractions() const
{
    return boost::make_iterator_range(d->improperTorsionInteractions.begin(),
                                      d->improperTorsionInteractions.end());
}

size_t Topology::improperTorsionInteractionCount() const
{
    return d->improperTorsionInteractions.size();
}

void Topology::addNonbondedInteraction(size_t i, size_t j)
{
    NonbondedInteraction interaction;
    interaction[0] = i;
    interaction[1] = j;
    d->nonbondedInteractions.push_back(interaction);
}

Topology::NonbondedInteractionRange Topology::nonbondedInteractions() const
{
    return boost::make_iterator_range(d->nonbondedInteractions.begin(),
                                      d->nonbondedInteractions.end());
}

size_t Topology::nonbondedInteractionCount() const
{
    return d->nonbondedInteractions.size();
}

/// Returns \c true if atoms \p i and \p j are in a one-four configuration.
bool Topology::isOneFour(size_t i, size_t j)
{
    foreach(const TorsionInteraction &torsion, d->torsionInteractions){
        if((torsion[0] == i && torsion[3] == j) ||
           (torsion[0] == j && torsion[3] == i)){
            return true;
        }
    }

    return false;
}

} // end chemkit namespace
