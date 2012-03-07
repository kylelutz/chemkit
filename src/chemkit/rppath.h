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

#ifndef CHEMKIT_RPPATH_H
#define CHEMKIT_RPPATH_H

#include "chemkit.h"

#include <set>
#include <limits>
#include <algorithm>

#include <boost/bind.hpp>

#include <Eigen/Core>

#include "atom.h"
#include "graph.h"
#include "foreach.h"
#include "fragment.h"
#include "molecule.h"

namespace chemkit {
namespace algorithm {
namespace detail {

// === PidMatrix =========================================================== //
// The PidMatrix class implements a path-included distance matrix.
template<typename T>
class PidMatrix
{
public:
    // construction and destruction
    PidMatrix(T size);
    ~PidMatrix();

    // paths
    std::vector<std::vector<T> >& paths(T i, T j);
    void addPaths(T i, T j, const std::vector<std::vector<T> > &paths);
    std::vector<std::vector<T> > splice(T i, T j, T k);

    // operators
    std::vector<std::vector<T> >& operator()(T i, T j);

private:
    T m_size;
    std::vector<std::vector<T> > *m_values;
};

// --- Construction and Destruction ---------------------------------------- //
template<typename T>
inline PidMatrix<T>::PidMatrix(T size)
{
    m_size = size;
    m_values = new std::vector<std::vector<T> >[size*size];
}

template<typename T>
inline PidMatrix<T>::~PidMatrix()
{
    delete [] m_values;
}

// --- Paths --------------------------------------------------------------- //
template<typename T>
inline std::vector<std::vector<T> >& PidMatrix<T>::paths(T i, T j)
{
    return m_values[i * m_size + j];
}

template<typename T>
inline void PidMatrix<T>::addPaths(T i, T j, const std::vector<std::vector<T> > &paths)
{
    std::vector<std::vector<T> > &current = m_values[i * m_size + j];
    current.insert(current.end(), paths.begin(), paths.end());
}

template<typename T>
inline std::vector<std::vector<T> >& PidMatrix<T>::operator()(T i, T j)
{
    return paths(i, j);
}

template<typename T>
inline std::vector<std::vector<T> > PidMatrix<T>::splice(T i, T j, T k)
{
    std::vector<std::vector<T> > splicedPaths;

    const std::vector<std::vector<T> > &ijPaths = paths(i, j);
    const std::vector<std::vector<T> > &jkPaths = paths(j, k);

    if(ijPaths.empty() && jkPaths.empty()){
        std::vector<T> path;
        path.push_back(j);
        splicedPaths.push_back(path);
    }
    else if(ijPaths.empty()){
        foreach(const std::vector<T> &jkPath, jkPaths){
            std::vector<T> path;
            path.push_back(j);
            path.insert(path.end(), jkPath.begin(), jkPath.end());
            splicedPaths.push_back(path);
        }
    }
    else if(jkPaths.empty()){
        foreach(const std::vector<T> &ijPath, ijPaths){
            std::vector<T> path = ijPath;
            path.push_back(j);
            splicedPaths.push_back(path);
        }
    }
    else{
        foreach(const std::vector<T> &ijPath, ijPaths){
            foreach(const std::vector<T> &jkPath, jkPaths){
                std::vector<T> path = ijPath;
                path.push_back(j);
                path.insert(path.end(), jkPath.begin(), jkPath.end());
                splicedPaths.push_back(path);
            }
        }
    }

    return splicedPaths;
}

// === RingCandidate ======================================================= //
template<typename T>
class RingCandidate
{
public:
    // construction and destruction
    RingCandidate(T size, T start, T end) : m_size(size), m_start(start), m_end(end) { }

    // properties
    T size() const { return m_size; }
    T start() const { return m_start; }
    T end() const { return m_end; }

    // static methods
    static bool compareSize(const RingCandidate &a, const RingCandidate &b) { return a.size() < b.size(); }

private:
    T m_size;
    T m_start;
    T m_end;
};

// === Sssr ================================================================ //
template<typename T>
class Sssr
{
public:
    // properties
    size_t size() const { return m_rings.size(); }
    bool isEmpty() const { return m_rings.empty(); }

    // rings
    const std::vector<std::vector<T> >& rings() const { return m_rings; }
    void append(const std::vector<T> &ring) { m_rings.push_back(ring); }

    // ring checks
    bool isValid(const std::vector<T> &ring) const;
    bool isUnique(const std::vector<T> &ring) const;

private:
    std::vector<std::vector<T> > m_rings;
};

// --- Ring Checks --------------------------------------------------------- //
template<typename T>
inline bool Sssr<T>::isValid(const std::vector<T> &ring) const
{
    // check for any duplicate verticies
    for(T i = 0; i < ring.size(); i++){
        for(T j = i + 1; j < ring.size(); j++){
            if(ring[i] == ring[j]){
                return false;
            }
        }
    }

    return true;
}

template<typename T>
inline bool Sssr<T>::isUnique(const std::vector<T> &path) const
{
    // must be unique if sssr is empty
    if(isEmpty()){
        return true;
    }

    // check if a ring with the same atoms is already in the sssr
    std::set<T> pathSet;
    pathSet.insert(path.begin(), path.end());

    foreach(const std::vector<T> &ring, m_rings){
        std::set<T> ringSet;
        ringSet.insert(ring.begin(), ring.end());

        std::vector<T> sortedRing(ring.begin(), ring.end());
        std::sort(sortedRing.begin(), sortedRing.end());

        std::set<T> intersection;
        std::set_intersection(pathSet.begin(), pathSet.end(),
                              ringSet.begin(), ringSet.end(),
                              std::inserter(intersection, intersection.begin()));

        if(intersection.size() == ring.size()){
            return false;
        }
    }

    // build set of bonds in the path
    std::set<std::pair<T, T> > pathBonds;
    for(T i = 0; i < path.size() - 1; i++){
        pathBonds.insert(std::make_pair(std::min(path[i], path[i+1]),
                                        std::max(path[i], path[i+1])));
    }

    pathBonds.insert(std::make_pair(std::min(path.front(), path.back()),
                                    std::max(path.front(), path.back())));

    // remove bonds from path bonds that are already in a smaller ring
    foreach(const std::vector<T> &ring, m_rings){
        if(ring.size() >= path.size()){
            continue;
        }

        for(T i = 0; i < ring.size() - 1; i++){
            pathBonds.erase(std::make_pair(std::min(ring[i], ring[i+1]),
                                           std::max(ring[i], ring[i+1])));
        }

        pathBonds.erase(std::make_pair(std::min(ring.front(), ring.back()),
                                       std::max(ring.front(), ring.back())));
    }

    // check if any other ring contains the same bonds
    foreach(const std::vector<T> &ring, m_rings){
        std::set<std::pair<T, T> > ringBonds;

        // add ring bonds
        for(T i = 0; i < ring.size() - 1; i++){
            ringBonds.insert(std::make_pair(std::min(ring[i], ring[i+1]),
                                            std::max(ring[i], ring[i+1])));
        }

        // add closure bond
        ringBonds.insert(std::make_pair(std::min(ring.front(), ring.back()),
                                        std::max(ring.front(), ring.back())));

        // check intersection
        std::set<std::pair<T, T> > intersection;
        std::set_intersection(pathBonds.begin(), pathBonds.end(),
                              ringBonds.begin(), ringBonds.end(),
                              std::inserter(intersection, intersection.begin()));

        if(intersection.size() == pathBonds.size()){
            return false;
        }
    }

    return true;
}

} // end detail namespace

// Returns the smallest set of smallest rings in a graph using the
// RP-Path algorithm.
//
// For a description of the algorithm see [Lee 2009].
template<typename T>
inline std::vector<std::vector<T> > rppath(const Graph<T> &graph)
{
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DistanceMatrix;
    using chemkit::algorithm::detail::PidMatrix;
    using chemkit::algorithm::detail::RingCandidate;
    using chemkit::algorithm::detail::Sssr;

    T n = graph.size();

    T ringCount = graph.edgeCount() - graph.vertexCount() + 1;
    if(ringCount == 0){
        return std::vector<std::vector<T> >();
    }

    // algorithm 1 - create the distance and pid matrices
    DistanceMatrix D(n, n);
    PidMatrix<T> P(n);
    PidMatrix<T> Pt(n);

    for(T i = 0; i < n; i++){
        for(T j = 0; j < n; j++){
            if(i == j){
                D(i, j) = 0;
            }
            else if(graph.isAdjacent(i, j)){
                D(i, j) = 1;
            }
            else{
                D(i, j) = std::numeric_limits<T>::max() / 2; // ~ infinity
            }
        }
    }

    for(T k = 0; k < n; k++){
        for(T i = 0; i < n; i++){
            for(T j = 0; j < n; j++){
                if(i == j || i == k || k == j){
                    continue;
                }

                if(D(i, j) > D(i, k) + D(k, j)){
                    if(D(i, j) == D(i, k) + D(k, j) + 1){
                        Pt(i, j) = P(i, j);
                    }
                    else{
                        Pt(i, j).clear();
                    }

                    D(i, j) = D(i, k) + D(k, j);
                    P(i, j) = P.splice(i, k, j);
                }
                else if(D(i, j) == D(i, k) + D(k, j)){
                    P.addPaths(i, j, P.splice(i, k, j));
                }
                else if(D(i, j) == D(i, k) + D(k, j) - 1){
                    Pt.addPaths(i, j, P.splice(i, k, j));
                }
            }
        }
    }

    // algorithm 2 - create the ring candidate set
    std::vector<RingCandidate<T> > candidates;
    for(T i = 0; i < n; i++){
        for(T j = i + 1; j < n; j++){
            if(P(i, j).size() == 1 && Pt(i, j).size() == 0){
                continue;
            }
            else{
                T size;

                if(P(i, j).size() > 1){
                    size = 2 * D(i, j);
                }
                else{
                    size = 2 * D(i, j) + 1;
                }

                if(size > 2){
                    candidates.push_back(RingCandidate<T>(size, i, j));
                }
            }
        }
    }

    // sort candidates
    std::sort(candidates.begin(), candidates.end(), RingCandidate<T>::compareSize);

    // algorithm 3 - find sssr from the ring candidate set
    Sssr<T> sssr;

    foreach(const RingCandidate<T> &candidate, candidates){
        // odd sized ring
        if(candidate.size() & 1){
            for(size_t i = 0; i < Pt(candidate.start(), candidate.end()).size(); i++){
                std::vector<T> ring;
                ring.push_back(candidate.start());
                std::vector<T> &path = Pt(candidate.start(), candidate.end())[i];
                ring.insert(ring.end(), path.begin(), path.end());
                ring.push_back(candidate.end());
                if(!P(candidate.end(), candidate.start()).empty()){
                    path = P(candidate.end(), candidate.start())[0];
                    ring.insert(ring.end(), path.begin(), path.end());
                }

                // check if ring is valid and unique
                if(sssr.isValid(ring) && sssr.isUnique(ring)){
                    sssr.append(ring);
                    break;
                }
            }
        }
        // even sized ring
        else{
            for(size_t i = 0; i < P(candidate.start(), candidate.end()).size() - 1; i++){
                std::vector<T> ring;
                ring.push_back(candidate.start());
                std::vector<T> &path = P(candidate.start(), candidate.end())[i];
                ring.insert(ring.end(), path.begin(), path.end());
                ring.push_back(candidate.end());
                path = P(candidate.end(), candidate.start())[i+1];
                ring.insert(ring.end(), path.begin(), path.end());

                // check if ring is valid and unique
                if(sssr.isValid(ring) && sssr.isUnique(ring)){
                    sssr.append(ring);
                    break;
                }
            }
        }

        if(sssr.size() == ringCount){
            break;
        }
    }

    return sssr.rings();
}

inline std::vector<std::vector<Atom *> > rppath(const Fragment *fragment)
{
    std::vector<Atom *> atoms = fragment->atoms();

    // remove any terminal atoms
    atoms.erase(std::remove_if(atoms.begin(), atoms.end(), boost::bind(&Atom::isTerminal, _1)), atoms.end());

    // create graph
    Graph<size_t> graph(atoms.size());

    for(size_t i = 0; i < atoms.size(); i++){
        for(size_t j = i + 1; j < atoms.size(); j++){
            if(atoms[i]->isBondedTo(atoms[j])){
                graph.addEdge(i, j);
            }
        }
    }

    // cyclize graph
    std::vector<size_t> originalIndicies;
    graph.cyclize(originalIndicies);

    // perceive rings
    std::vector<std::vector<size_t> > sssr = rppath(graph);

    // convert from lists of indices to lists of atoms
    std::vector<std::vector<Atom *> > rings;

    foreach(const std::vector<size_t> &cycle, sssr){
        std::vector<Atom *> ring(cycle.size());

        for(size_t i = 0; i < cycle.size(); i++){
            ring[i] = atoms[originalIndicies[cycle[i]]];
        }

        rings.push_back(ring);
    }

    return rings;
}

inline std::vector<std::vector<Atom *> > rppath(const Molecule *molecule)
{
    std::vector<std::vector<Atom *> > rings;

    foreach(const Fragment *fragment, molecule->fragments()){
        std::vector<std::vector<Atom *> > fragmentRings = rppath(fragment);

        rings.insert(rings.end(), fragmentRings.begin(), fragmentRings.end());
    }

    return rings;
}

} // end algorithm namespace
} // end chemkit namespace

#endif // CHEMKIT_ALGORITHM_RPPATH_H
