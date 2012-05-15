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

#ifndef CHEMKIT_MOLECULEGRAPHTRAITS_H
#define CHEMKIT_MOLECULEGRAPHTRAITS_H

#include <utility>

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include "atom.h"
#include "bond.h"
#include "molecule.h"

namespace boost {

// === Molecule Graph Edge Iterator ======================================== //
class chemkit_molecule_graph_edge_iterator :
    public iterator_facade<chemkit_molecule_graph_edge_iterator,
                           remove_reference<
                               const std::pair<const chemkit::Bond*, bool> >::type,
                           multi_pass_input_iterator_tag,
                           remove_reference<
                               const std::pair<const chemkit::Bond*, bool> >::type>
{
public:
    typedef chemkit::Atom::BondRange::const_iterator BaseIterator;

    chemkit_molecule_graph_edge_iterator()
    {
    }

    chemkit_molecule_graph_edge_iterator(const BaseIterator &base,
                                         const chemkit::Atom *source)
        : m_iterator(base),
          m_source(source)
    {
    }

    chemkit_molecule_graph_edge_iterator(
        const chemkit_molecule_graph_edge_iterator &other)
        : m_iterator(other.m_iterator),
          m_source(other.m_source)
    {
    }

    chemkit_molecule_graph_edge_iterator& operator=(
        const chemkit_molecule_graph_edge_iterator &other)
    {
        m_iterator = other.m_iterator;
        m_source = other.m_source;

        return *this;
    }

private:
    bool equal(const chemkit_molecule_graph_edge_iterator &other) const
    {
        return m_iterator == other.m_iterator;
    }

    void advance(size_t distance)
    {
        std::advance(m_iterator, distance);
    }

    void increment()
    {
        m_iterator++;
    }

    void decrement()
    {
        m_iterator--;
    }

    const std::pair<const chemkit::Bond *, bool> dereference() const
    {
        return std::make_pair(*m_iterator, m_source == (*m_iterator)->atom1());
    }

    friend class iterator_core_access;

private:
    BaseIterator m_iterator;
    const chemkit::Atom *m_source;
};

// === Molecule Graph Traits =============================================== //
struct chemkit_molecule_graph_traversal_category :
        public virtual incidence_graph_tag,
        public virtual bidirectional_graph_tag,
        public virtual adjacency_graph_tag,
        public virtual vertex_list_graph_tag,
        public virtual edge_list_graph_tag
{
};

template<>
struct graph_traits<chemkit::Molecule>
{
    typedef const chemkit::Atom* vertex_descriptor;
    typedef std::pair<const chemkit::Bond*, bool> edge_descriptor;
    typedef chemkit::Atom::NeighborRange::const_iterator adjacency_iterator;
    typedef chemkit_molecule_graph_edge_iterator out_edge_iterator;
    typedef chemkit_molecule_graph_edge_iterator in_edge_iterator;
    typedef chemkit::Molecule::AtomRange::const_iterator vertex_iterator;
    typedef chemkit_molecule_graph_edge_iterator edge_iterator;
    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category;
    typedef chemkit_molecule_graph_traversal_category traversal_category;
    typedef size_t vertices_size_type;
    typedef size_t edges_size_type;
    typedef size_t degree_size_type;

    static vertex_descriptor null_vertex() { return 0; }
};

inline const chemkit::Atom*
source(graph_traits<chemkit::Molecule>::edge_descriptor edge,
       const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    const chemkit::Bond *bond = edge.first;
    bool forward = edge.second;

    return forward ? bond->atom1() : bond->atom2();
}

inline const chemkit::Atom*
target(graph_traits<chemkit::Molecule>::edge_descriptor edge,
       const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    const chemkit::Bond *bond = edge.first;
    bool forward = edge.second;

    return forward ? bond->atom2() : bond->atom1();
}

inline std::pair<graph_traits<chemkit::Molecule>::out_edge_iterator,
                 graph_traits<chemkit::Molecule>::out_edge_iterator>
out_edges(graph_traits<chemkit::Molecule>::vertex_descriptor vertex,
          const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    const chemkit::Atom::BondRange bonds = vertex->bonds();

    return std::make_pair(chemkit_molecule_graph_edge_iterator(bonds.begin(), vertex),
                          chemkit_molecule_graph_edge_iterator(bonds.end(), vertex));
}

inline std::pair<graph_traits<chemkit::Molecule>::in_edge_iterator,
                 graph_traits<chemkit::Molecule>::in_edge_iterator>
in_edges(graph_traits<chemkit::Molecule>::vertex_descriptor vertex,
         const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    const chemkit::Atom::BondRange bonds = vertex->bonds();

    return std::make_pair(chemkit_molecule_graph_edge_iterator(bonds.begin(), vertex),
                          chemkit_molecule_graph_edge_iterator(bonds.end(), vertex));
}

inline graph_traits<chemkit::Molecule>::degree_size_type
out_degree(graph_traits<chemkit::Molecule>::vertex_descriptor vertex,
           const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    return vertex->neighborCount();
}

inline graph_traits<chemkit::Molecule>::degree_size_type
in_degree(graph_traits<chemkit::Molecule>::vertex_descriptor vertex,
          const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    return vertex->neighborCount();
}

inline std::pair<graph_traits<chemkit::Molecule>::vertex_iterator,
                 graph_traits<chemkit::Molecule>::vertex_iterator>
vertices(const chemkit::Molecule &graph)
{
    const chemkit::Molecule::AtomRange atoms = graph.atoms();

    return std::make_pair(atoms.begin(), atoms.end());
}

inline graph_traits<chemkit::Molecule>::vertices_size_type
num_vertices(const chemkit::Molecule &graph)
{
    return graph.atomCount();
}

inline std::pair<graph_traits<chemkit::Molecule>::adjacency_iterator,
                 graph_traits<chemkit::Molecule>::adjacency_iterator>
adjacent_vertices(graph_traits<chemkit::Molecule>::vertex_descriptor vertex,
                  const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    const chemkit::Atom::NeighborRange neighbors = vertex->neighbors();

    return std::make_pair(neighbors.begin(), neighbors.end());
}

inline std::pair<graph_traits<chemkit::Molecule>::edge_iterator,
                 graph_traits<chemkit::Molecule>::edge_iterator>
edges(const chemkit::Molecule &graph)
{
    const chemkit::Molecule::BondRange bonds = graph.bonds();

    return std::make_pair(chemkit_molecule_graph_edge_iterator(bonds.begin(), 0),
                          chemkit_molecule_graph_edge_iterator(bonds.end(), 0));
}

inline graph_traits<chemkit::Molecule>::edges_size_type
num_edges(const chemkit::Molecule &graph)
{
    return graph.bondCount();
}

// === Molecule Graph Property Maps ======================================== //
struct chemkit_molecule_graph_vertex_index_map
{
public:
    typedef size_t value_type;
    typedef size_t reference;
    typedef graph_traits<chemkit::Molecule>::vertex_descriptor key_type;
    typedef readable_property_map_tag category;
};

inline size_t
get(chemkit_molecule_graph_vertex_index_map,
    graph_traits<chemkit::Molecule>::vertex_descriptor vertex)
{
    return vertex->index();
}

inline chemkit_molecule_graph_vertex_index_map
get(vertex_index_t, const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    return chemkit_molecule_graph_vertex_index_map();
}

template<>
struct property_map<chemkit::Molecule, vertex_index_t>
{
    typedef chemkit_molecule_graph_vertex_index_map type;
    typedef chemkit_molecule_graph_vertex_index_map const_type;
};

struct chemkit_molecule_graph_edge_index_map
{
public:
    typedef size_t value_type;
    typedef size_t reference;
    typedef graph_traits<chemkit::Molecule>::edge_descriptor key_type;
    typedef readable_property_map_tag category;
};

inline size_t
get(chemkit_molecule_graph_edge_index_map,
    graph_traits<chemkit::Molecule>::edge_descriptor edge)
{
    return edge.first->index();
}

inline chemkit_molecule_graph_edge_index_map
get(edge_index_t, const chemkit::Molecule &graph)
{
    CHEMKIT_UNUSED(graph);

    return chemkit_molecule_graph_edge_index_map();
}

template<>
struct property_map<chemkit::Molecule, edge_index_t>
{
    typedef chemkit_molecule_graph_edge_index_map type;
    typedef chemkit_molecule_graph_edge_index_map const_type;
};

} // end boost namespace

#endif // CHEMKIT_MOLECULEGRAPHTRAITS_H
