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

#ifndef CHEMKIT_GRAPH_INLINE_H
#define CHEMKIT_GRAPH_INLINE_H

#include "graph.h"

#include <cassert>
#include <algorithm>

#include "foreach.h"

namespace chemkit {

// === Graph =============================================================== //
/// \class Graph graph.h chemkit/graph.h
/// \ingroup chemkit
/// \internal
/// \brief The Graph class represents a graph.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graph with \p size vertices.
template<typename T>
inline Graph<T>::Graph(T size)
    : m_adjacencyList(size)
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the number of vertices in the graph to \p size.
template<typename T>
inline void Graph<T>::resize(T size)
{
    m_adjacencyList.resize(size);
}

/// Returns the number of vertices in the graph.
template<typename T>
inline T Graph<T>::size() const
{
    return vertexCount();
}

/// Returns \c true if the graph contains no vertices.
template<typename T>
inline bool Graph<T>::isEmpty() const
{
    return m_adjacencyList.empty();
}

// --- Structure ----------------------------------------------------------- //
/// Adds a new vertex to the graph and returns its index.
template<typename T>
inline T Graph<T>::addVertex()
{
    resize(size() + 1);

    return size() - 1;
}

/// Removes \p vertex from the graph.
template<typename T>
inline void Graph<T>::removeVertex(T vertex)
{
    m_adjacencyList.erase(m_adjacencyList.begin() + vertex);
}

/// Returns the number of vertices in the graph.
template<typename T>
inline T Graph<T>::vertexCount() const
{
    return m_adjacencyList.size();
}

/// Adds an edge between vertices \p a and \p b.
template<typename T>
inline void Graph<T>::addEdge(T a, T b)
{
    assert(!isAdjacent(a, b));

    m_adjacencyList[a].push_back(b);
    m_adjacencyList[b].push_back(a);
}

/// Removes the edge between vertices \p a and \p b.
template<typename T>
inline void Graph<T>::removeEdge(T a, T b)
{
    assert(isAdjacent(a, b));

    m_adjacencyList[a].erase(std::find(m_adjacencyList[a].begin(), m_adjacencyList[a].end(), b));
    m_adjacencyList[b].erase(std::find(m_adjacencyList[b].begin(), m_adjacencyList[b].end(), a));
}

/// Returns the number of edges in the graph.
template<typename T>
inline T Graph<T>::edgeCount() const
{
    T count = 0;

    for(size_t i = 0; i < m_adjacencyList.size(); i++){
        count += m_adjacencyList[i].size();
    }

    return count / 2;
}

/// Returns \c true if vertex \p a is adjacent to vertex \p b.
template<typename T>
inline bool Graph<T>::isAdjacent(T a, T b) const
{
    const std::vector<T> &neighborsA = neighbors(a);

    return std::find(neighborsA.begin(), neighborsA.end(), b) != neighborsA.end();
}

// --- Algorithms ---------------------------------------------------------- //
/// Swap vertices \p a and \p b.
template<typename T>
inline void Graph<T>::swap(T a, T b)
{
    // update a's neighbors
    foreach(T neighbor, m_adjacencyList[a]){
        size_t index = std::distance(m_adjacencyList[neighbor].begin(),
                                     std::find(m_adjacencyList[neighbor].begin(),
                                               m_adjacencyList[neighbor].end(),
                                               a));
        m_adjacencyList[neighbor][index] = b;
    }

    // update b's neighbors
    foreach(T neighbor, m_adjacencyList[b]){
        size_t index = std::distance(m_adjacencyList[neighbor].begin(),
                                     std::find(m_adjacencyList[neighbor].begin(),
                                               m_adjacencyList[neighbor].end(),
                                               b));
        m_adjacencyList[neighbor][index] = a;
    }

    // swap the adjacency lists for each vertex
    std::swap(m_adjacencyList[a], m_adjacencyList[b]);
}

/// Returns the neighbors of vertex.
template<typename T>
inline const std::vector<T>& Graph<T>::neighbors(T vertex) const
{
    return m_adjacencyList[vertex];
}

/// Remove all terminal vertices from the graph.
template<typename T>
inline void Graph<T>::cyclize(std::vector<T> &originalIndices)
{
    // remove all terminal edges
    bool done = false;
    while(!done){
        done = true;

        for(size_t i = 0; i < m_adjacencyList.size(); i++){
            if(m_adjacencyList[i].size() == 1){
                T neighbor = m_adjacencyList[i][0];
                removeEdge(i, neighbor);
                done = false;
            }
        }
    }

    originalIndices.resize(m_adjacencyList.size());
    for(size_t i = 0; i < m_adjacencyList.size(); i++){
        originalIndices[i] = i;
    }

    for(size_t i = 0; i < m_adjacencyList.size(); i++){
        // find lone vertex
        if(m_adjacencyList[i].empty()){
            // find next non-terminal vertex
            for(size_t j = i + 1; j < m_adjacencyList.size(); j++){
                if(!m_adjacencyList[j].empty()){
                    // swap the terminal and non-terminal vertices
                    swap(i, j);
                    originalIndices[i] = j;
                    break;
                }
            }
        }
    }

    // remove the lone vertices
    for(size_t i = 0; i < m_adjacencyList.size(); i++){
        if(m_adjacencyList[i].empty()){
            m_adjacencyList.resize(i);
            originalIndices.resize(i);
            break;
        }
    }
}

} // end chemkit namespace

#endif // CHEMKIT_GRAPH_INLINE_H
