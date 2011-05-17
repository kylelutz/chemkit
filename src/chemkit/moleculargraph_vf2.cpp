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

// This file implements the VF2 graph isomorphism algorithm. See
// [Cordella 2005].

#include "moleculargraph.h"

#include "foreach.h"

namespace chemkit {

namespace {

// The SharedState class holds four vectors containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
class SharedState
{
    public:
        SharedState(unsigned int sourceSize, unsigned int targetSize);

        std::vector<int> sourceMapping;
        std::vector<int> targetMapping;
        std::vector<unsigned int> sourceTerminalSet;
        std::vector<unsigned int> targetTerminalSet;
};

SharedState::SharedState(unsigned int sourceSize, unsigned int targetSize)
    : sourceMapping(sourceSize, -1),
      targetMapping(targetSize, -1),
      sourceTerminalSet(sourceSize, 0),
      targetTerminalSet(targetSize, 0)
{
}

// The State class represents a single state in the isomorphism detection
// algorithm. Every state uses and modifies the same SharedState object.
class State
{
    public:
        State(const MolecularGraph *source, const MolecularGraph *target);
        State(const State *state);
        ~State();

        unsigned int size() const { return m_size; }
        const MolecularGraph* source() const { return m_source; }
        const MolecularGraph* target() const { return m_target; }
        const Atom* sourceAtom(int index) { return m_source->atom(index); }
        const Atom* targetAtom(int index) { return m_target->atom(index); }
        std::map<Atom *, Atom *> mapping() const;
        bool succeeded() const;
        void addPair(const std::pair<int, int> &candidate);
        std::pair<int, int> nextCandidate(const std::pair<int, int> &lastCandidate);
        bool isFeasible(const std::pair<int, int> &candidate);
        void backTrack();

    private:
        unsigned int m_size;
        unsigned int m_sourceTerminalSize;
        unsigned int m_targetTerminalSize;
        const MolecularGraph *m_source;
        const MolecularGraph *m_target;
        std::pair<int, int> m_lastAddition;
        SharedState *m_sharedState;
        bool m_ownSharedState;
};

State::State(const MolecularGraph *source, const MolecularGraph *target)
    : m_size(0),
      m_sourceTerminalSize(0),
      m_targetTerminalSize(0),
      m_source(source),
      m_target(target),
      m_lastAddition(-1, -1),
      m_sharedState(new SharedState(source->size(), target->size())),
      m_ownSharedState(true)
{
}

State::State(const State *state)
    : m_size(state->m_size),
      m_sourceTerminalSize(state->m_sourceTerminalSize),
      m_targetTerminalSize(state->m_targetTerminalSize),
      m_source(state->m_source),
      m_target(state->m_target),
      m_lastAddition(-1, -1),
      m_sharedState(state->m_sharedState),
      m_ownSharedState(false)
{
}

State::~State()
{
    if(m_ownSharedState)
        delete m_sharedState;
}

// Returns true if the state contains an isomorphism.
bool State::succeeded() const
{
    return m_size == m_source->size();
}

// Returns the current isomorphism for the state as a std::map.
std::map<Atom *, Atom *> State::mapping() const
{
    std::map<Atom *, Atom *> mapping;

    for(unsigned int i = 0; i < m_size; i++){
        mapping[m_source->atom(i)] = m_target->atom(m_sharedState->sourceMapping[i]);
    }

    return mapping;
}

// Returns the next candidate pair (sourceAtom, targetAtom) to be added to the
// state. The candidate should be checked for feasibility and then added using
// the addPair() method.
std::pair<int, int> State::nextCandidate(const std::pair<int, int> &lastCandidate)
{
    int lastSourceAtom = lastCandidate.first;
    int lastTargetAtom = lastCandidate.second;

    int sourceSize = m_source->size();
    int targetSize = m_target->size();

    if(lastSourceAtom == -1)
        lastSourceAtom = 0;

    if(lastTargetAtom == -1)
        lastTargetAtom = 0;
    else
        lastTargetAtom++;

    if(m_sourceTerminalSize > m_size && m_targetTerminalSize > m_size){
        while(lastSourceAtom < sourceSize &&
              (m_sharedState->sourceMapping[lastSourceAtom] != -1 ||
               m_sharedState->sourceTerminalSet[lastSourceAtom] == 0)){
            lastSourceAtom++;
            lastTargetAtom = 0;
        }
    }
    else{
        while(lastSourceAtom < sourceSize &&
              m_sharedState->sourceMapping[lastSourceAtom] != -1){
            lastSourceAtom++;
            lastTargetAtom = 0;
        }
    }

    if(m_sourceTerminalSize > m_size && m_targetTerminalSize > m_size){
        while(lastTargetAtom < targetSize &&
              (m_sharedState->targetMapping[lastTargetAtom] != -1 ||
               m_sharedState->targetTerminalSet[lastTargetAtom] == 0)){
            lastTargetAtom++;
        }
    }
    else{
        while(lastTargetAtom < targetSize &&
              m_sharedState->targetMapping[lastTargetAtom] != -1){
            lastTargetAtom++;
        }
    }

    if(lastSourceAtom < sourceSize && lastTargetAtom < targetSize){
        return std::make_pair(lastSourceAtom, lastTargetAtom);
    }

    return std::make_pair(-1, -1);
}

// Adds the candidate pair (sourceAtom, targetAtom) to the state. The candidate
// pair must be feasible to add it to the state.
void State::addPair(const std::pair<int, int> &candidate)
{
    m_size++;
    m_lastAddition = candidate;

    int sourceAtom = candidate.first;
    int targetAtom = candidate.second;

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

    foreach(int neighbor, m_source->neighbors(sourceAtom)){
        if(!m_sharedState->sourceTerminalSet[neighbor]){
            m_sharedState->sourceTerminalSet[neighbor] = m_size;
            m_sourceTerminalSize++;
        }
    }

    foreach(int neighbor, m_target->neighbors(targetAtom)){
        if(!m_sharedState->targetTerminalSet[neighbor]){
            m_sharedState->targetTerminalSet[neighbor] = m_size;
            m_targetTerminalSize++;
        }
    }
}

// Restores the shared state to how it was before adding the last candidate
// pair. Assumes addPair() has been called on the state only once.
void State::backTrack()
{
    int addedSourceAtom = m_lastAddition.first;

    if(m_sharedState->sourceTerminalSet[addedSourceAtom] == m_size){
        m_sharedState->sourceTerminalSet[addedSourceAtom] = 0;
    }

    foreach(int neighbor, m_source->neighbors(addedSourceAtom)){
        if(m_sharedState->sourceTerminalSet[neighbor] == m_size){
            m_sharedState->sourceTerminalSet[neighbor] = 0;
        }
    }

    int addedTargetAtom = m_lastAddition.second;

    if(m_sharedState->targetTerminalSet[addedTargetAtom] == m_size){
        m_sharedState->targetTerminalSet[addedTargetAtom] = 0;
    }

    foreach(int neighbor, m_target->neighbors(addedTargetAtom)){
        if(m_sharedState->targetTerminalSet[neighbor] == m_size){
            m_sharedState->targetTerminalSet[neighbor] = 0;
        }
    }

    m_sharedState->sourceMapping[addedSourceAtom] = -1;
    m_sharedState->targetMapping[addedTargetAtom] = -1;
    m_size--;
    m_lastAddition = std::make_pair(-1, -1);
}

bool State::isFeasible(const std::pair<int, int> &candidate)
{
    int sourceAtom = candidate.first;
    int targetAtom = candidate.second;

    int sourceAtomLabel = m_source->atomLabel(sourceAtom);
    int targetAtomLabel = m_target->atomLabel(targetAtom);

    if(sourceAtomLabel != targetAtomLabel){
        return false;
    }

    int sourceTerminalNeighborCount = 0;
    int targetTerminalNeighborCount = 0;
    int sourceNewNeighborCount = 0;
    int targetNewNeighborCount = 0;

    foreach(int neighbor, m_source->neighbors(sourceAtom)){
        int sourceBond = m_source->bond(sourceAtom, neighbor);
        int sourceBondLabel = m_source->bondLabel(sourceBond);

        if(m_sharedState->sourceMapping[neighbor] != -1){
            int targetNeighbor = m_sharedState->sourceMapping[neighbor];

            if(!m_target->isAdjacent(targetAtom, targetNeighbor))
                return false;

            int targetBond = m_target->bond(targetAtom, targetNeighbor);
            int targetBondLabel = m_target->bondLabel(targetBond);

            if(sourceBondLabel != targetBondLabel)
                return false;
        }
        else{
            if(m_sharedState->sourceTerminalSet[neighbor])
                sourceTerminalNeighborCount++;
            else
                sourceNewNeighborCount++;
        }
    }

    foreach(int neighbor, m_target->neighbors(targetAtom)){
        if(m_sharedState->targetMapping[neighbor] != -1){
            //int sourceNeighbor = m_sharedState->targetMapping[neighbor];

            //if(!m_source->adjacent(sourceAtom, sourceNeighbor)){
            //    return false;
            //}
        }
        else{
            if(m_sharedState->targetTerminalSet[neighbor])
                targetTerminalNeighborCount++;
            else
                targetNewNeighborCount++;
        }
    }

    return (sourceTerminalNeighborCount <= targetTerminalNeighborCount) &&
           (sourceNewNeighborCount <= targetNewNeighborCount);
}

bool match(State *state, std::map<Atom *, Atom *> &mapping)
{
    if(state->succeeded()){
        mapping = state->mapping();
        return true;
    }

    std::pair<int, int> lastCanidate(-1, -1);

    bool found = false;
    while(!found){
        std::pair<int, int> candidate = state->nextCandidate(lastCanidate);

        if(candidate.first == -1)
            return false;

        lastCanidate = candidate;

        if(state->isFeasible(candidate)){
            State nextState(state);
            nextState.addPair(candidate);
            found = match(&nextState, mapping);
            nextState.backTrack();
        }
    }

    return found;
}

} // end anonymous namespace

// The isomorphism_vf2() method returns an isomorphism between two molecular
// graphs using the VF2 algorithm. This can be used for finding both
// graph-graph isomorphisms and graph-subgraph isomorphisms. In the latter case
// graph 'a' is the subgraph, implying a->size() < b->size(). In the case that
// no isomorphism is found an empty mapping is returned.
std::map<Atom *, Atom *> MolecularGraph::isomorphism_vf2(const MolecularGraph *a, const MolecularGraph *b)
{
    State state(a, b);
    std::map<Atom *, Atom *> mapping;
    match(&state, mapping);
    return mapping;
}

} // end chemkit namespace
