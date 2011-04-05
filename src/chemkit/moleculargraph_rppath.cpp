/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

// This file implements the RP-Path ring perception algorithm. See
// [Lee 2009].

#include "moleculargraph.h"

#include <set>
#include <algorithm>

#include "ring.h"
#include "foreach.h"

namespace chemkit {

namespace {

// === DistanceMatrix ====================================================== //
class DistanceMatrix
{
    public:
        // construction and destruction
        DistanceMatrix(int size);
        ~DistanceMatrix();

        // operators
        int operator()(int i, int j) const;
        int& operator()(int i, int j);

    private:
        int m_size;
        int *m_values;
};

DistanceMatrix::DistanceMatrix(int size)
{
    m_size = size;
    m_values = new int[size*size];
    memset(m_values, 0, size*size*sizeof(int));
}

DistanceMatrix::~DistanceMatrix()
{
    delete [] m_values;
}

int DistanceMatrix::operator()(int i, int j) const
{
    return m_values[i * m_size + j];
}

int& DistanceMatrix::operator()(int i, int j)
{
    return m_values[i * m_size + j];
}

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

    std::vector<std::vector<int> > ijPaths = paths(i, j);
    std::vector<std::vector<int> > jkPaths = paths(j, k);

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

// --- Functions ----------------------------------------------------------- //
bool isValid(const std::vector<int> &path)
{
    // check for any duplicate atoms
    for(unsigned int i = 0; i < path.size(); i++){
        for(unsigned int j = i + 1; j < path.size(); j++){
            if(path[i] == path[j]){
                return false;
            }
        }
    }

    return true;
}

bool isUnique(const std::vector<int> &path, const std::vector<std::vector<int> > &sssr)
{
    // check if a ring with the same atoms is already in the sssr
    std::vector<int> sortedPath(path.begin(), path.end());
    std::sort(sortedPath.begin(), sortedPath.end());

    foreach(const std::vector<int> &ring, sssr){
        if(path.size() != ring.size()){
            continue;
        }

        std::vector<int> sortedRing(ring.begin(), ring.end());
        std::sort(sortedRing.begin(), sortedRing.end());

        bool different = false;
        for(unsigned int i = 0; i < path.size(); i++){
            if(sortedPath[i] != sortedRing[i]){
                different = true;
                break;
            }
        }

        if(!different){
            return false;
        }
    }

    // count number of unique bonds
    std::set<std::pair<int, int> > ringBonds;
    foreach(const std::vector<int> &ring, sssr){

        // add ring bonds
        for(unsigned int i = 0; i < ring.size()-1; i++){
            ringBonds.insert(std::make_pair(qMin(ring[i], ring[i+1]),
                                            qMax(ring[i], ring[i+1])));
        }

        // add closure bond
        ringBonds.insert(std::make_pair(qMin(ring.front(), ring.back()),
                                        qMax(ring.front(), ring.back())));
    }

    int uniqueBondCount = 0;
    for(unsigned int i = 0; i < path.size()-1; i++){
        if(ringBonds.find(std::make_pair(qMin(path[i], path[i+1]),
                                         qMax(path[i], path[i+1]))) == ringBonds.end()){
            uniqueBondCount++;
        }
    }

    // count closure bond
    if(ringBonds.find(std::make_pair(qMin(path.front(), path.back()),
                                     qMax(path.front(), path.back()))) == ringBonds.end()){
        uniqueBondCount++;
    }

    return uniqueBondCount > 0;
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
    DistanceMatrix D(n);
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
                D(i, j) = INT_MAX/2; // ~ infinity
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
    std::vector<std::vector<int> > sssr;
    sssr.reserve(ringCount);

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
                if(isValid(ring) && isUnique(ring, sssr)){
                    sssr.push_back(ring);
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
                if(isValid(ring) && isUnique(ring, sssr)){
                    sssr.push_back(ring);
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
    foreach(const std::vector<int> &ring, sssr){
        std::vector<Atom *> atoms;

        foreach(int atomIndex, ring){
            atoms.push_back(graph->atom(atomIndex));
        }

        rings.push_back(new Ring(atoms));
    }

    return rings;
}

} // end chemkit namespace
