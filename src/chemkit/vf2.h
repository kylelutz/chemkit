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

#ifndef CHEMKIT_VF2_H
#define CHEMKIT_VF2_H

#include "chemkit.h"

#include <map>
#include <vector>

#include "graph.h"

namespace chemkit {
namespace algorithm {
namespace detail {

// The SharedState class holds four vectors containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
template<typename T>
class SharedState
{
public:
    typedef T SizeType;
    enum { NullIndex = SizeType(-1) }; // represents an invalid vertex index

    SharedState(T sourceSize, T targetSize);

    std::vector<T> sourceMapping;
    std::vector<T> targetMapping;
    std::vector<T> sourceTerminalSet;
    std::vector<T> targetTerminalSet;
};

template<typename T>
inline SharedState<T>::SharedState(T sourceSize, T targetSize)
    : sourceMapping(sourceSize, NullIndex),
      targetMapping(targetSize, NullIndex),
      sourceTerminalSet(sourceSize, 0),
      targetTerminalSet(targetSize, 0)
{
}

// The State class represents a single state in the isomorphism detection
// algorithm. Every state uses and modifies the same SharedState object.
template<typename T, typename VertexComparator, typename EdgeComparator>
class State
{
public:
    typedef T SizeType;
    enum { NullIndex = SizeType(-1) }; // represents an invalid vertex index

    State(const Graph<T> &source, const Graph<T> &target, VertexComparator compareVertices, EdgeComparator compareEdges);
    State(const State *state);
    ~State();

    SizeType size() const { return m_size; }
    const Graph<T>& source() const { return m_source; }
    const Graph<T>& target() const { return m_target; }
    std::map<T, T> mapping() const;
    bool succeeded() const;
    void addPair(const std::pair<T, T> &candidate);
    std::pair<T, T> nextCandidate(const std::pair<T, T> &lastCandidate);
    bool isFeasible(const std::pair<T, T> &candidate);
    void backTrack();

    static const std::pair<T, T> nullCandidate() { return std::make_pair(SizeType(-1), SizeType(-1)); }

private:
    SizeType m_size;
    SizeType m_sourceTerminalSize;
    SizeType m_targetTerminalSize;
    const Graph<T> &m_source;
    const Graph<T> &m_target;
    std::pair<T, T> m_lastAddition;
    SharedState<T> *m_sharedState;
    bool m_ownSharedState;
    VertexComparator m_compareVertices;
    EdgeComparator m_compareEdges;
};

template<typename T, typename VertexComparator, typename EdgeComparator>
inline State<T, VertexComparator, EdgeComparator>::State(const Graph<T> &source,
                                                         const Graph<T> &target,
                                                         VertexComparator compareVertices,
                                                         EdgeComparator compareEdges)
    : m_size(0),
      m_sourceTerminalSize(0),
      m_targetTerminalSize(0),
      m_source(source),
      m_target(target),
      m_lastAddition(NullIndex, NullIndex),
      m_sharedState(new SharedState<T>(source.size(), target.size())),
      m_ownSharedState(true),
      m_compareVertices(compareVertices),
      m_compareEdges(compareEdges)
{
}

template<typename T, typename VertexComparator, typename EdgeComparator>
inline State<T, VertexComparator, EdgeComparator>::State(const State *state)
    : m_size(state->m_size),
      m_sourceTerminalSize(state->m_sourceTerminalSize),
      m_targetTerminalSize(state->m_targetTerminalSize),
      m_source(state->m_source),
      m_target(state->m_target),
      m_lastAddition(NullIndex, NullIndex),
      m_sharedState(state->m_sharedState),
      m_ownSharedState(false),
      m_compareVertices(state->m_compareVertices),
      m_compareEdges(state->m_compareEdges)
{
}

template<typename T, typename VertexComparator, typename EdgeComparator>
inline State<T, VertexComparator, EdgeComparator>::~State()
{
    if(m_ownSharedState)
        delete m_sharedState;
}

// Returns true if the state contains an isomorphism.
template<typename T, typename VertexComparator, typename EdgeComparator>
inline bool State<T, VertexComparator, EdgeComparator>::succeeded() const
{
    return m_size == m_source.size();
}

// Returns the current isomorphism for the state as a std::map.
template<typename T, typename VertexComparator, typename EdgeComparator>
inline std::map<T, T> State<T, VertexComparator, EdgeComparator>::mapping() const
{
    std::map<T, T> mapping;

    for(SizeType i = 0; i < m_size; i++){
        mapping[i] = m_sharedState->sourceMapping[i];
    }

    return mapping;
}

// Returns the next candidate pair (sourceAtom, targetAtom) to be added to the
// state. The candidate should be checked for feasibility and then added using
// the addPair() method.
template<typename T, typename VertexComparator, typename EdgeComparator>
inline std::pair<T, T> State<T, VertexComparator, EdgeComparator>::nextCandidate(const std::pair<T, T> &lastCandidate)
{
    T lastSourceAtom = lastCandidate.first;
    T lastTargetAtom = lastCandidate.second;

    SizeType sourceSize = m_source.size();
    SizeType targetSize = m_target.size();

    if(lastSourceAtom == NullIndex)
        lastSourceAtom = 0;

    if(lastTargetAtom == NullIndex)
        lastTargetAtom = 0;
    else
        lastTargetAtom++;

    if(m_sourceTerminalSize > m_size && m_targetTerminalSize > m_size){
        while(lastSourceAtom < sourceSize &&
              (m_sharedState->sourceMapping[lastSourceAtom] != NullIndex ||
               m_sharedState->sourceTerminalSet[lastSourceAtom] == 0)){
            lastSourceAtom++;
            lastTargetAtom = 0;
        }
    }
    else{
        while(lastSourceAtom < sourceSize &&
              m_sharedState->sourceMapping[lastSourceAtom] != NullIndex){
            lastSourceAtom++;
            lastTargetAtom = 0;
        }
    }

    if(m_sourceTerminalSize > m_size && m_targetTerminalSize > m_size){
        while(lastTargetAtom < targetSize &&
              (m_sharedState->targetMapping[lastTargetAtom] != NullIndex||
               m_sharedState->targetTerminalSet[lastTargetAtom] == 0)){
            lastTargetAtom++;
        }
    }
    else{
        while(lastTargetAtom < targetSize &&
              m_sharedState->targetMapping[lastTargetAtom] != NullIndex){
            lastTargetAtom++;
        }
    }

    if(lastSourceAtom < sourceSize && lastTargetAtom < targetSize){
        return std::make_pair(lastSourceAtom, lastTargetAtom);
    }

    return nullCandidate();
}

// Adds the candidate pair (sourceAtom, targetAtom) to the state. The candidate
// pair must be feasible to add it to the state.
template<typename T, typename VertexComparator, typename EdgeComparator>
inline void State<T, VertexComparator, EdgeComparator>::addPair(const std::pair<T, T> &candidate)
{
    m_size++;
    m_lastAddition = candidate;

    T sourceAtom = candidate.first;
    T targetAtom = candidate.second;

    if(!m_sharedState->sourceTerminalSet[sourceAtom]){
        m_sharedState->sourceTerminalSet[sourceAtom] = m_size;
        //m_sourceTerminalSize++;
    }

    if(!m_sharedState->targetTerminalSet[targetAtom]){
        m_sharedState->targetTerminalSet[targetAtom] = m_size;
        //m_targetTerminalSize++;
    }

    m_sharedState->sourceMapping[sourceAtom] = targetAtom;
    m_sharedState->targetMapping[targetAtom] = sourceAtom;

    foreach(T neighbor, m_source.neighbors(sourceAtom)){
        if(!m_sharedState->sourceTerminalSet[neighbor]){
            m_sharedState->sourceTerminalSet[neighbor] = m_size;
            m_sourceTerminalSize++;
        }
    }

    foreach(T neighbor, m_target.neighbors(targetAtom)){
        if(!m_sharedState->targetTerminalSet[neighbor]){
            m_sharedState->targetTerminalSet[neighbor] = m_size;
            m_targetTerminalSize++;
        }
    }
}

// Restores the shared state to how it was before adding the last candidate
// pair. Assumes addPair() has been called on the state only once.
template<typename T, typename VertexComparator, typename EdgeComparator>
inline void State<T, VertexComparator, EdgeComparator>::backTrack()
{
    T addedSourceAtom = m_lastAddition.first;

    if(m_sharedState->sourceTerminalSet[addedSourceAtom] == m_size){
        m_sharedState->sourceTerminalSet[addedSourceAtom] = 0;
    }

    foreach(T neighbor, m_source.neighbors(addedSourceAtom)){
        if(m_sharedState->sourceTerminalSet[neighbor] == m_size){
            m_sharedState->sourceTerminalSet[neighbor] = 0;
        }
    }

    T addedTargetAtom = m_lastAddition.second;

    if(m_sharedState->targetTerminalSet[addedTargetAtom] == m_size){
        m_sharedState->targetTerminalSet[addedTargetAtom] = 0;
    }

    foreach(T neighbor, m_target.neighbors(addedTargetAtom)){
        if(m_sharedState->targetTerminalSet[neighbor] == m_size){
            m_sharedState->targetTerminalSet[neighbor] = 0;
        }
    }

    m_sharedState->sourceMapping[addedSourceAtom] = NullIndex;
    m_sharedState->targetMapping[addedTargetAtom] = NullIndex;
    m_size--;
    m_lastAddition = nullCandidate();
}

template<typename T, typename VertexComparator, typename EdgeComparator>
inline bool State<T, VertexComparator, EdgeComparator>::isFeasible(const std::pair<T, T> &candidate)
{
    T sourceAtom = candidate.first;
    T targetAtom = candidate.second;

    if(!m_compareVertices(sourceAtom, targetAtom)){
        return false;
    }

    SizeType sourceTerminalNeighborCount = 0;
    SizeType targetTerminalNeighborCount = 0;
    SizeType sourceNewNeighborCount = 0;
    SizeType targetNewNeighborCount = 0;

    foreach(T neighbor, m_source.neighbors(sourceAtom)){
        if(m_sharedState->sourceMapping[neighbor] != NullIndex){
            T targetNeighbor = m_sharedState->sourceMapping[neighbor];

            if(!m_target.isAdjacent(targetAtom, targetNeighbor)){
                return false;
            }

            if(!m_compareEdges(sourceAtom, neighbor, targetAtom, targetNeighbor)){
                return false;
            }
        }
        else{
            if(m_sharedState->sourceTerminalSet[neighbor]){
                sourceTerminalNeighborCount++;
            }
            else{
                sourceNewNeighborCount++;
            }
        }
    }

    foreach(T neighbor, m_target.neighbors(targetAtom)){
        if(m_sharedState->targetMapping[neighbor] == NullIndex){
            if(m_sharedState->targetTerminalSet[neighbor]){
                targetTerminalNeighborCount++;
            }
            else{
                targetNewNeighborCount++;
            }
        }
    }

    return (sourceTerminalNeighborCount <= targetTerminalNeighborCount) &&
           (sourceNewNeighborCount <= targetNewNeighborCount);
}

template<typename T, typename VertexComparator, typename EdgeComparator>
inline bool match(State<T, VertexComparator, EdgeComparator> *state, std::map<T, T> &mapping)
{
    if(state->succeeded()){
        mapping = state->mapping();
        return true;
    }

    std::pair<T, T> lastCandidate = state->nullCandidate();

    bool found = false;
    while(!found){
        std::pair<T, T> candidate = state->nextCandidate(lastCandidate);

        if(candidate == state->nullCandidate()){
            return false;
        }

        lastCandidate = candidate;

        if(state->isFeasible(candidate)){
            State<T, VertexComparator, EdgeComparator> nextState(state);
            nextState.addPair(candidate);
            found = match(&nextState, mapping);
            nextState.backTrack();
        }
    }

    return found;
}

} // end detail namespace

template<typename T, typename VertexComparator, typename EdgeComparator>
std::map<T, T> vf2(const Graph<T> &a,
                   const Graph<T> &b,
                   VertexComparator vertexComparator,
                   EdgeComparator edgeComparator)
{
    using detail::State;
    using detail::match;

    // create initial empty state
    State<T, VertexComparator, EdgeComparator> state(a,
                                                     b,
                                                     vertexComparator,
                                                     edgeComparator);

    // create empty mapping
    std::map<T, T> mapping;

    // run vf2 match algorithm
    match(&state, mapping);

    // return mapping
    return mapping;
}

} // end algorithm namespace
} // end chemkit namespace

#endif // CHEMKIT_VF2_H
