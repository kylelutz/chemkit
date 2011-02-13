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

// This file implements the VF2 graph isomorphism algorithm. See
// [Cordella 2005].

#include "moleculargraph.h"

#include "molecule.h"
#include "atommapping.h"

namespace chemkit {

namespace {

// The SharedState class holds four vectors containing the mapping between
// the two graphs and the terminal sets. It is shared between all the states
// in each isomorphism test.
class SharedState
{
    public:
        SharedState(int sourceSize, int targetSize);

        QVector<int> sourceMapping;
        QVector<int> targetMapping;
        QVector<int> sourceTerminalSet;
        QVector<int> targetTerminalSet;
};

SharedState::SharedState(int sourceSize, int targetSize)
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

        int size() const { return m_size; }
        const MolecularGraph* source() const { return m_source; }
        const MolecularGraph* target() const { return m_target; }
        const Atom* sourceAtom(int index) { return m_source->atom(index); }
        const Atom* targetAtom(int index) { return m_target->atom(index); }
        AtomMapping mapping() const;
        bool succeeded() const;
        void addPair(const QPair<int, int> &candidate);
        QPair<int, int> nextCandidate(const QPair<int, int> &lastCandidate);
        bool isFeasible(const QPair<int, int> &candidate);
        void backTrack();

    private:
        int m_size;
        int m_sourceTerminalSize;
        int m_targetTerminalSize;
        const MolecularGraph *m_source;
        const MolecularGraph *m_target;
        QPair<int, int> m_lastAddition;
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

// Returns the current isomorphism for the state in an AtomMapping object.
AtomMapping State::mapping() const
{
    AtomMapping mapping(m_source->molecule(), m_target->molecule());

    for(int i = 0; i < m_size; i++){
        mapping.add(m_source->atom(i), m_target->atom(m_sharedState->sourceMapping[i]));
    }

    return mapping;
}

// Returns the next candidate pair (sourceAtom, targetAtom) to be added to the
// state. The candidate should be checked for feasibility and then added using
// the addPair() method.
QPair<int, int> State::nextCandidate(const QPair<int, int> &lastCandidate)
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
        return qMakePair(lastSourceAtom, lastTargetAtom);
    }

    return qMakePair(-1, -1);
}

// Adds the candidate pair (sourceAtom, targetAtom) to the state. The candidate
// pair must be feasible to add it to the state.
void State::addPair(const QPair<int, int> &candidate)
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
    m_lastAddition = qMakePair(-1, -1);
}

bool State::isFeasible(const QPair<int, int> &candidate)
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

            if(!m_target->adjacent(targetAtom, targetNeighbor))
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

bool match(State *state, AtomMapping &mapping)
{
    if(state->succeeded()){
        mapping = state->mapping();
        return true;
    }

    QPair<int, int> lastCanidate(-1, -1);

    bool found = false;
    while(!found){
        QPair<int, int> candidate = state->nextCandidate(lastCanidate);

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
AtomMapping MolecularGraph::isomorphism_vf2(const MolecularGraph *a, const MolecularGraph *b)
{
    State state(a, b);
    AtomMapping mapping(a->molecule(), b->molecule());
    match(&state, mapping);
    return mapping;
}

} // end chemkit namespace
