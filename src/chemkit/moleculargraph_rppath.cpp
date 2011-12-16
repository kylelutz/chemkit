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

// This file implements the RP-Path ring perception algorithm. See
// [Lee 2009].

#include "moleculargraph.h"

#include <set>
#include <limits>
#include <algorithm>

#include <Eigen/Core>

#include "ring.h"
#include "foreach.h"

namespace chemkit {

namespace {

typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> DistanceMatrix;

// === PidMatrix =========================================================== //
// The PidMatrix class implements a path-included distance matrix.
class PidMatrix
{
public:
    // construction and destruction
    PidMatrix(int size);
    ~PidMatrix();

    // paths
    std::vector<std::vector<int> >& paths(int i, int j);
    void addPaths(int i, int j, const std::vector<std::vector<int> > &paths);
    std::vector<std::vector<int> > splice(int i, int j, int k);

    // operators
    std::vector<std::vector<int> >& operator()(int i, int j);

private:
    int m_size;
    std::vector<std::vector<int> > *m_values;
};

// --- Construction and Destruction ---------------------------------------- //
PidMatrix::PidMatrix(int size)
{
    m_size = size;
    m_values = new std::vector<std::vector<int> >[size*size];
}

PidMatrix::~PidMatrix()
{
    delete [] m_values;
}

// --- Paths --------------------------------------------------------------- //
std::vector<std::vector<int> >& PidMatrix::paths(int i, int j)
{
    return m_values[i * m_size + j];
}

void PidMatrix::addPaths(int i, int j, const std::vector<std::vector<int> > &paths)
{
    std::vector<std::vector<int> > &current = m_values[i * m_size + j];
    current.insert(current.end(), paths.begin(), paths.end());
}

std::vector<std::vector<int> >& PidMatrix::operator()(int i, int j)
{
    return paths(i, j);
}

std::vector<std::vector<int> > PidMatrix::splice(int i, int j, int k)
{
    std::vector<std::vector<int> > splicedPaths;

    const std::vector<std::vector<int> > &ijPaths = paths(i, j);
    const std::vector<std::vector<int> > &jkPaths = paths(j, k);

    if(ijPaths.empty() && jkPaths.empty()){
        std::vector<int> path;
        path.push_back(j);
        splicedPaths.push_back(path);
    }
    else if(ijPaths.empty()){
        foreach(const std::vector<int> &jkPath, jkPaths){
            std::vector<int> path;
            path.push_back(j);
            path.insert(path.end(), jkPath.begin(), jkPath.end());
            splicedPaths.push_back(path);
        }
    }
    else if(jkPaths.empty()){
        foreach(const std::vector<int> &ijPath, ijPaths){
            std::vector<int> path = ijPath;
            path.push_back(j);
            splicedPaths.push_back(path);
        }
    }
    else{
        foreach(const std::vector<int> &ijPath, ijPaths){
            foreach(const std::vector<int> &jkPath, jkPaths){
                std::vector<int> path = ijPath;
                path.push_back(j);
                path.insert(path.end(), jkPath.begin(), jkPath.end());
                splicedPaths.push_back(path);
            }
        }
    }

    return splicedPaths;
}

// === RingCandidate ======================================================= //
class RingCandidate
{
public:
    // construction and destruction
    RingCandidate(int size, int start, int end);

    // properties
    int size() const;
    int start() const;
    int end() const;

    // static methods
    static bool compareSize(const RingCandidate &a, const RingCandidate &b);

private:
    int m_size;
    int m_start;
    int m_end;
};

// --- Construction and Destruction ---------------------------------------- //
RingCandidate::RingCandidate(int size, int start, int end)
{
    m_size = size;
    m_start = start;
    m_end = end;
}

// --- Properties ---------------------------------------------------------- //
int RingCandidate::size() const
{
    return m_size;
}

int RingCandidate::start() const
{
    return m_start;
}

int RingCandidate::end() const
{
    return m_end;
}

// --- Static Methods ------------------------------------------------------ //
bool RingCandidate::compareSize(const RingCandidate &a, const RingCandidate &b)
{
    return a.size() < b.size();
}

// === Sssr ================================================================ //
class Sssr
{
public:
    // construction and destruction
    Sssr();
    ~Sssr();

    // properties
    unsigned int size() const;
    bool isEmpty() const;

    // rings
    const std::vector<std::vector<int> >& rings() const;
    void append(const std::vector<int> &ring);
    bool isValid(const std::vector<int> &ring) const;
    bool isUnique(const std::vector<int> &ring) const;

private:
    std::vector<std::vector<int> > m_rings;
};

// --- Construction and Destruction ---------------------------------------- //
Sssr::Sssr()
{
}

Sssr::~Sssr()
{
}

// --- Properties ---------------------------------------------------------- //
unsigned int Sssr::size() const
{
    return m_rings.size();
}

bool Sssr::isEmpty() const
{
    return m_rings.empty();
}

// --- Rings --------------------------------------------------------------- //
const std::vector<std::vector<int> >& Sssr::rings() const
{
    return m_rings;
}

void Sssr::append(const std::vector<int> &ring)
{
    m_rings.push_back(ring);
}

bool Sssr::isValid(const std::vector<int> &ring) const
{
    // check for any duplicate atoms
    for(unsigned int i = 0; i < ring.size(); i++){
        for(unsigned int j = i + 1; j < ring.size(); j++){
            if(ring[i] == ring[j]){
                return false;
            }
        }
    }

    return true;
}

bool Sssr::isUnique(const std::vector<int> &path) const
{
    // must be unique if sssr is empty
    if(isEmpty()){
        return true;
    }

    // check if a ring with the same atoms is already in the sssr
    std::set<int> pathSet;
    pathSet.insert(path.begin(), path.end());

    foreach(const std::vector<int> &ring, m_rings){
        std::set<int> ringSet;
        ringSet.insert(ring.begin(), ring.end());

        std::vector<int> sortedRing(ring.begin(), ring.end());
        std::sort(sortedRing.begin(), sortedRing.end());

        std::set<int> intersection;
        std::set_intersection(pathSet.begin(), pathSet.end(),
                              ringSet.begin(), ringSet.end(),
                              std::inserter(intersection, intersection.begin()));

        if(intersection.size() == ring.size()){
            return false;
        }
    }

    // build set of bonds in the path
    std::set<std::pair<int, int> > pathBonds;
    for(unsigned int i = 0; i < path.size()-1; i++){
        pathBonds.insert(std::make_pair(std::min(path[i], path[i+1]),
                                         std::max(path[i], path[i+1])));
    }

    pathBonds.insert(std::make_pair(std::min(path.front(), path.back()),
                                    std::max(path.front(), path.back())));

    // remove bonds from path bonds that are already in a smaller ring
    foreach(const std::vector<int> &ring, m_rings){
        if(ring.size() >= path.size()){
            continue;
        }

        for(unsigned int i = 0; i < ring.size() - 1; i++){
            pathBonds.erase(std::make_pair(std::min(ring[i], ring[i+1]),
                                           std::max(ring[i], ring[i+1])));
        }

        pathBonds.erase(std::make_pair(std::min(ring.front(), ring.back()),
                                       std::max(ring.front(), ring.back())));
    }

    // check if any other ring contains the same bonds
    foreach(const std::vector<int> &ring, m_rings){
        std::set<std::pair<int, int> > ringBonds;

        // add ring bonds
        for(unsigned int i = 0; i < ring.size()-1; i++){
            ringBonds.insert(std::make_pair(std::min(ring[i], ring[i+1]),
                                            std::max(ring[i], ring[i+1])));
        }

        // add closure bond
        ringBonds.insert(std::make_pair(std::min(ring.front(), ring.back()),
                                        std::max(ring.front(), ring.back())));

        // check intersection
        std::set<std::pair<int, int> > intersection;
        std::set_intersection(pathBonds.begin(), pathBonds.end(),
                              ringBonds.begin(), ringBonds.end(),
                              std::inserter(intersection, intersection.begin()));

        if(intersection.size() == pathBonds.size()){
            return false;
        }
    }

    return true;
}

} // end anonymous namespace

// The sssr_rpPath() method returns the smallest set of smallest rings in a
// molecular graph using the RP-Path algorithm. The graph is expected to
// contain a single fragment and to have all terminal nodes removed (i.e. all
// verticies should have degree >= 2).
std::vector<Ring *> MolecularGraph::sssr_rpPath(const MolecularGraph *graph)
{
    unsigned int n = graph->size();

    unsigned int ringCount = graph->bondCount() - graph->atomCount() + 1;
    if(ringCount == 0){
        return std::vector<Ring *>();
    }

    // algorithm 1 - create the distance and pid matrices
    DistanceMatrix D(n, n);
    PidMatrix P(n);
    PidMatrix Pt(n);

    for(unsigned int i = 0; i < n; i++){
        for(unsigned int j = 0; j < n; j++){
            if(i == j){
                D(i, j) = 0;
            }
            else if(graph->isAdjacent(i, j)){
                D(i, j) = 1;
            }
            else{
                D(i, j) = std::numeric_limits<unsigned int>::max() / 2; // ~ infinity
            }
        }
    }

    for(unsigned int k = 0; k < n; k++){
        for(unsigned int i = 0; i < n; i++){
            for(unsigned int j = 0; j < n; j++){
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
    std::vector<RingCandidate> candidates;
    for(unsigned int i = 0; i < n; i++){
        for(unsigned int j = i + 1; j < n; j++){
            if(P(i, j).size() == 1 && Pt(i, j).size() == 0){
                continue;
            }
            else{
                int size;

                if(P(i, j).size() > 1){
                    size = 2 * D(i, j);
                }
                else{
                    size = 2 * D(i, j) + 1;
                }

                if(size > 2){
                    candidates.push_back(RingCandidate(size, i, j));
                }
            }
        }
    }

    // sort candidates
    std::sort(candidates.begin(), candidates.end(), RingCandidate::compareSize);

    // algorithm 3 - find sssr from the ring candidate set
    Sssr sssr;

    foreach(const RingCandidate &candidate, candidates){

        // odd sized ring
        if(candidate.size() & 1){
            for(unsigned int i = 0; i < Pt(candidate.start(), candidate.end()).size(); i++){
                std::vector<int> ring;
                ring.push_back(candidate.start());
                std::vector<int> &path = Pt(candidate.start(), candidate.end())[i];
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
            for(unsigned int i = 0; i < P(candidate.start(), candidate.end()).size()-1; i++){
                std::vector<int> ring;
                ring.push_back(candidate.start());
                std::vector<int> &path = P(candidate.start(), candidate.end())[i];
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

    // build list of rings
    std::vector<Ring *> rings;
    foreach(const std::vector<int> &ring, sssr.rings()){
        std::vector<Atom *> atoms;

        foreach(int atomIndex, ring){
            atoms.push_back(graph->atom(atomIndex));
        }

        rings.push_back(new Ring(atoms));
    }

    return rings;
}

} // end chemkit namespace
