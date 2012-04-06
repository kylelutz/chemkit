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

#include "substructurequery.h"

#include <boost/make_shared.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>

#include "vf2.h"
#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "foreach.h"
#include "molecule.h"

namespace chemkit {

namespace {

struct AtomComparator
{
    AtomComparator(const std::vector<Atom *> &sourceAtoms, const std::vector<Atom *> &targetAtoms)
        : m_sourceAtoms(sourceAtoms),
          m_targetAtoms(targetAtoms)
    {
    }

    AtomComparator(const AtomComparator &other)
        : m_sourceAtoms(other.m_sourceAtoms),
          m_targetAtoms(other.m_targetAtoms)
    {
    }

    bool operator()(size_t a, size_t b) const
    {
        return m_sourceAtoms[a]->atomicNumber() == m_targetAtoms[b]->atomicNumber();
    }

    const std::vector<Atom *> &m_sourceAtoms;
    const std::vector<Atom *> &m_targetAtoms;
};

struct BondComparator
{
    BondComparator(const std::vector<Atom *> &sourceAtoms, const std::vector<Atom *> &targetAtoms, int flags)
        : m_sourceAtoms(sourceAtoms),
          m_targetAtoms(targetAtoms),
          m_flags(flags)
    {
    }

    BondComparator(const BondComparator &other)
        : m_sourceAtoms(other.m_sourceAtoms),
          m_targetAtoms(other.m_targetAtoms),
          m_flags(other.m_flags)
    {
    }

    bool operator()(size_t a1, size_t a2, size_t b1, size_t b2) const
    {
        const Bond *bondA = m_sourceAtoms[a1]->bondTo(m_sourceAtoms[a2]);
        const Bond *bondB = m_targetAtoms[b1]->bondTo(m_targetAtoms[b2]);

        if(!bondA || !bondB){
            return false;
        }

        if(m_flags & SubstructureQuery::CompareAromaticity){
            return (bondA->order() == bondB->order()) ||
                   (bondA->isAromatic() && bondB->isAromatic());
        }
        else{
            return bondA->order() == bondB->order();
        }
    }

    const std::vector<Atom *> &m_sourceAtoms;
    const std::vector<Atom *> &m_targetAtoms;
    int m_flags;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> AdjacencyListGraph;

struct AdjacencyListGraphVertexComparator
{
    AdjacencyListGraphVertexComparator(const std::vector<Atom *> &sourceAtoms,
                                       const std::vector<Atom *> &targetAtoms)
        : m_sourceAtoms(sourceAtoms),
          m_targetAtoms(targetAtoms)
    {
    }

    AdjacencyListGraphVertexComparator(const AdjacencyListGraphVertexComparator &other)
        : m_sourceAtoms(other.m_sourceAtoms),
          m_targetAtoms(other.m_targetAtoms)
    {
    }

    bool operator()(const boost::graph_traits<AdjacencyListGraph>::vertex_descriptor &vertexA,
                    const boost::graph_traits<AdjacencyListGraph>::vertex_descriptor &vertexB)
    {
        return m_sourceAtoms[vertexA]->atomicNumber() == m_targetAtoms[vertexB]->atomicNumber();
    }

    const std::vector<Atom *> &m_sourceAtoms;
    const std::vector<Atom *> &m_targetAtoms;
};

struct AdjacencyListGraphEdgeComparator
{
    AdjacencyListGraphEdgeComparator(const AdjacencyListGraph &source,
                                     const AdjacencyListGraph &target,
                                     const std::vector<Atom *> &sourceAtoms,
                                     const std::vector<Atom *> &targetAtoms,
                                     int flags)
        : m_source(source),
          m_target(target),
          m_sourceAtoms(sourceAtoms),
          m_targetAtoms(targetAtoms),
          m_flags(flags)
    {
    }

    AdjacencyListGraphEdgeComparator(const AdjacencyListGraphEdgeComparator &other)
        : m_source(other.m_source),
          m_target(other.m_target),
          m_sourceAtoms(other.m_sourceAtoms),
          m_targetAtoms(other.m_targetAtoms),
          m_flags(other.m_flags)
    {
    }

    bool operator()(const boost::graph_traits<AdjacencyListGraph>::edge_descriptor &edgeA,
                    const boost::graph_traits<AdjacencyListGraph>::edge_descriptor &edgeB)
    {
        size_t a1 = boost::source(edgeA, m_source);
        size_t a2 = boost::target(edgeA, m_source);
        size_t b1 = boost::source(edgeB, m_target);
        size_t b2 = boost::target(edgeB, m_target);

        const Bond *bondA = m_sourceAtoms[a1]->bondTo(m_sourceAtoms[a2]);
        const Bond *bondB = m_targetAtoms[b1]->bondTo(m_targetAtoms[b2]);

        if(!bondA || !bondB){
            return false;
        }

        if(m_flags & SubstructureQuery::CompareAromaticity){
            return (bondA->order() == bondB->order()) ||
                   (bondA->isAromatic() && bondB->isAromatic());
        }
        else{
            return bondA->order() == bondB->order();
        }
    }

    const AdjacencyListGraph &m_source;
    const AdjacencyListGraph &m_target;
    const std::vector<Atom *> &m_sourceAtoms;
    const std::vector<Atom *> &m_targetAtoms;
    int m_flags;
};

struct McgregorCommonSubgraphsCallback
{
    McgregorCommonSubgraphsCallback(const AdjacencyListGraph &source,
                                    const AdjacencyListGraph &target,
                                    std::map<size_t, size_t> &mapping)
        : m_source(source),
          m_target(target),
          m_mapping(mapping)
    {
    }

    McgregorCommonSubgraphsCallback(const McgregorCommonSubgraphsCallback &other)
        : m_source(other.m_source),
          m_target(other.m_target),
          m_mapping(other.m_mapping)
    {
    }

    std::map<size_t, size_t> mapping() const
    {
        return m_mapping;
    }

    template<typename MapType>
    bool operator()(MapType map1To2,
                    MapType map2To1,
                    typename boost::graph_traits<AdjacencyListGraph>::vertices_size_type size)
    {
        CHEMKIT_UNUSED(map2To1);
        CHEMKIT_UNUSED(size);

        if(m_mapping.empty()){
            BGL_FORALL_VERTICES_T(vertex1, m_source, AdjacencyListGraph){
                if(boost::get(map1To2, vertex1) != boost::graph_traits<AdjacencyListGraph>::null_vertex()) {
                    m_mapping[vertex1] = boost::get(map1To2, vertex1);
                }
            }
        }

        return false;
    }

    const AdjacencyListGraph &m_source;
    const AdjacencyListGraph &m_target;
    std::map<size_t, size_t> &m_mapping;
};

} // end anonymous namespace

// === SubstructureQueryPrivate ============================================ //
class SubstructureQueryPrivate
{
public:
    boost::shared_ptr<Molecule> molecule;
    int flags;
};

// === SubstructureQuery =================================================== //
/// \class SubstructureQuery substructurequery.h chemkit/substructurequery.h
/// \ingroup chemkit
/// \brief The SubstructureQuery class represents a substructure query.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new substructure query.
SubstructureQuery::SubstructureQuery()
  : d(new SubstructureQueryPrivate)
{
    d->flags = 0;
}

/// Creates a new substructure query with \p molecule as the
/// substructure to query for.
SubstructureQuery::SubstructureQuery(const boost::shared_ptr<Molecule> &molecule)
    : d(new SubstructureQueryPrivate)
{
    d->molecule = molecule;
    d->flags = 0;
}

/// Creates a new substructure query with \p formula in \p format as
/// the substructure to query for.
SubstructureQuery::SubstructureQuery(const std::string &formula, const std::string &format)
    : d(new SubstructureQueryPrivate)
{
    d->molecule = boost::make_shared<Molecule>(formula, format);
    d->flags = 0;
}

/// Destroys the substructure query object.
SubstructureQuery::~SubstructureQuery()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the substructure molecule to \p molecule.
void SubstructureQuery::setMolecule(const boost::shared_ptr<Molecule> &molecule)
{
    d->molecule = molecule;
}

/// Sets the substructure molecule to \p formula with \p format.
void SubstructureQuery::setMolecule(const std::string &formula, const std::string &format)
{
    setMolecule(boost::make_shared<Molecule>(formula, format));
}

/// Returns the substructure molecule.
boost::shared_ptr<Molecule> SubstructureQuery::molecule() const
{
    return d->molecule;
}

/// Sets the query flags to \p flags.
void SubstructureQuery::setFlags(int flags)
{
    d->flags = flags;
}

/// Returns the query flags.
int SubstructureQuery::flags() const
{
    return d->flags;
}

// --- Queries ------------------------------------------------------------- //
/// Returns \c true if the substructure molecule matches \p molecule.
///
/// For example, to find if a molecule contains a carbonyl group (C=O):
/// \code
/// bool containsCarbonylGroup(const Molecule *molecule)
/// {
///      boost::shared_ptr<Molecule> carbonyl(new Molecule);
///      Atom *C1 = carbonyl->addAtom("C");
///      Atom *O2 = carbonyl->addAtom("O");
///      carbonyl->addBond(C1, O2, Bond::Double);
///
///      SubstructureQuery query(carbonyl);
///
///      return query.matches(molecule);
/// }
/// \endcode
bool SubstructureQuery::matches(const Molecule *molecule) const
{
    if(!d->molecule){
        return false;
    }

    if(d->molecule->isEmpty()){
        return true;
    }

    return !mapping(molecule).empty();
}

/// Returns a mapping (also known as an isomorphism) between the
/// atoms in the substructure molecule and the atoms in \p molecule.
std::map<Atom *, Atom *> SubstructureQuery::mapping(const Molecule *molecule) const
{
    Graph<size_t> source;
    Graph<size_t> target;

    std::vector<Atom *> sourceAtoms;
    std::vector<Atom *> targetAtoms;

    if(d->flags & CompareHydrogens){
        source.resize(d->molecule->size());
        sourceAtoms = std::vector<Atom *>(d->molecule->atoms().begin(), d->molecule->atoms().end());

        if(!(d->flags & CompareAtomsOnly)){
            foreach(const Bond *bond, d->molecule->bonds()){
                source.addEdge(bond->atom1()->index(), bond->atom2()->index());
            }
        }

        target.resize(molecule->size());
        targetAtoms = std::vector<Atom *>(molecule->atoms().begin(), molecule->atoms().end());

        if(!(d->flags & CompareAtomsOnly)){
            foreach(const Bond *bond, molecule->bonds()){
                target.addEdge(bond->atom1()->index(), bond->atom2()->index());
            }
        }
    }
    else{
        foreach(Atom *atom, d->molecule->atoms()){
            if(!atom->isTerminalHydrogen()){
                sourceAtoms.push_back(atom);
            }
        }

        foreach(Atom *atom, molecule->atoms()){
            if(!atom->isTerminalHydrogen()){
                targetAtoms.push_back(atom);
            }
        }

        source.resize(sourceAtoms.size());
        target.resize(targetAtoms.size());

        if(!(d->flags & CompareAtomsOnly)){
            for(size_t i = 0; i < sourceAtoms.size(); i++){
                for(size_t j = i + 1; j < sourceAtoms.size(); j++){
                    if(sourceAtoms[i]->isBondedTo(sourceAtoms[j])){
                        source.addEdge(i, j);
                    }
                }
            }

            for(size_t i = 0; i < targetAtoms.size(); i++){
                for(size_t j = i + 1; j < targetAtoms.size(); j++){
                    if(targetAtoms[i]->isBondedTo(targetAtoms[j])){
                        target.addEdge(i, j);
                    }
                }
            }
        }
    }

    AtomComparator atomComparator(sourceAtoms, targetAtoms);
    BondComparator bondComparator(sourceAtoms, targetAtoms, d->flags);

    // run vf2 isomorphism algorithm
    std::map<size_t, size_t> mapping = chemkit::algorithm::vf2(source,
                                                               target,
                                                               atomComparator,
                                                               bondComparator);

    // check for exact match
    if(d->flags & CompareExact && mapping.size() != source.size()){
        return std::map<Atom *, Atom *>();
    }

    // convert index mapping to an atom mapping
    std::map<Atom *, Atom *> atomMapping;

    for(std::map<size_t, size_t>::iterator i = mapping.begin(); i != mapping.end(); ++i){
        atomMapping[sourceAtoms[i->first]] = targetAtoms[i->second];
    }

    return atomMapping;
}

/// Returns the maximum mapping (also known as maximum common
/// substructure or MCS) between the query molecule and \p molecule.
std::map<Atom *, Atom *> SubstructureQuery::maximumMapping(const Molecule *molecule) const
{
    AdjacencyListGraph source;
    AdjacencyListGraph target;

    std::vector<Atom *> sourceAtoms;
    std::vector<Atom *> targetAtoms;

    if(d->flags & CompareHydrogens){
        for(size_t i = 0; i < d->molecule->size(); i++){
            boost::add_vertex(source);
        }

        sourceAtoms = std::vector<Atom *>(d->molecule->atoms().begin(), d->molecule->atoms().end());

        if(!(d->flags & CompareAtomsOnly)){
            foreach(const Bond *bond, d->molecule->bonds()){
                boost::add_edge(bond->atom1()->index(), bond->atom2()->index(), source);
            }
        }

        for(size_t i = 0; i < molecule->size(); i++){
            boost::add_vertex(target);
        }

        targetAtoms = std::vector<Atom *>(molecule->atoms().begin(), molecule->atoms().end());

        if(!(d->flags & CompareAtomsOnly)){
            foreach(const Bond *bond, molecule->bonds()){
                boost::add_edge(bond->atom1()->index(), bond->atom2()->index(), target);
            }
        }
    }
    else{
        foreach(Atom *atom, d->molecule->atoms()){
            if(!atom->isTerminalHydrogen()){
                sourceAtoms.push_back(atom);
                boost::add_vertex(source);
            }
        }

        foreach(Atom *atom, molecule->atoms()){
            if(!atom->isTerminalHydrogen()){
                targetAtoms.push_back(atom);
                boost::add_vertex(target);
            }
        }

        if(!(d->flags & CompareAtomsOnly)){
            for(size_t i = 0; i < sourceAtoms.size(); i++){
                for(size_t j = i + 1; j < sourceAtoms.size(); j++){
                    if(sourceAtoms[i]->isBondedTo(sourceAtoms[j])){
                        boost::add_edge(i, j, source);
                    }
                }
            }

            for(size_t i = 0; i < targetAtoms.size(); i++){
                for(size_t j = i + 1; j < targetAtoms.size(); j++){
                    if(targetAtoms[i]->isBondedTo(targetAtoms[j])){
                        boost::add_edge(i, j, target);
                    }
                }
            }
        }
    }

    AdjacencyListGraphVertexComparator vertexComparator(sourceAtoms, targetAtoms);
    AdjacencyListGraphEdgeComparator edgeComparator(source, target, sourceAtoms, targetAtoms, d->flags);

    std::map<size_t, size_t> mapping;
    McgregorCommonSubgraphsCallback callback(source, target, mapping);

    // search for connected subgraphs if the query molecule
    // consists only of a single connected component
    bool onlyConnectedSubgraphs = !d->molecule->isFragmented();

    boost::mcgregor_common_subgraphs_maximum_unique(source,
                                                    target,
                                                    boost::get(boost::vertex_index, source),
                                                    boost::get(boost::vertex_index, target),
                                                    edgeComparator,
                                                    vertexComparator,
                                                    onlyConnectedSubgraphs,
                                                    callback);

    // convert index mapping to an atom mapping
    std::map<Atom *, Atom *> atomMapping;

    for(std::map<size_t, size_t>::iterator i = mapping.begin(); i != mapping.end(); ++i){
        atomMapping[sourceAtoms[i->first]] = targetAtoms[i->second];
    }

    return atomMapping;
}

/// Returns a vector containing each molecule in \p molecules that
/// matches the substructure molecule.
std::vector<Molecule *> SubstructureQuery::filter(const std::vector<Molecule *> &molecules) const
{
    std::vector<Molecule *> matchingMolecules;

    foreach(Molecule *molecule, molecules){
        if(matches(molecule)){
            matchingMolecules.push_back(molecule);
        }
    }

    return matchingMolecules;
}

/// Searches the the molecule for an occurrence of the substructure
/// molecule in \p molecule and returns it if found. If not found an
/// empty moiety is returned.
///
/// For example, to find an amide group (NC=O) in a molecule:
/// \code
/// boost::shared_ptr<Molecule> amide(new Molecule);
/// Atom *C1 = amide->addAtom("C");
/// Atom *N2 = amide->addAtom("N");
/// Atom *O3 = amide->addAtom("O");
/// amide->addBond(C1, N2, Bond::Single);
/// amide->addBond(C1, O3, Bond::Double);
///
/// SubstructureQuery query(amide);
///
/// Moiety amideGroup = query.find(molecule);
/// \endcode
Moiety SubstructureQuery::find(const Molecule *molecule) const
{
    std::map<Atom *, Atom *> mapping = this->mapping(molecule);

    // no mapping found, return empty moiety
    if(mapping.empty()){
        return Moiety();
    }

    std::vector<Atom *> atoms;

    foreach(Atom *atom, d->molecule->atoms()){
        atoms.push_back(mapping[atom]);
    }

    return Moiety(atoms);
}

} // end chemkit namespace
