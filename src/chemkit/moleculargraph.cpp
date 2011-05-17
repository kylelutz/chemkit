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

#include "moleculargraph.h"

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "foreach.h"
#include "molecule.h"

#include <algorithm>

namespace chemkit {

// === MolecularGraph ====================================================== //
/// \class MolecularGraph moleculargraph.h chemkit/moleculargraph.h
/// \ingroup chemkit
/// \internal
/// \brief The MolecularGraph class represents a molecular graph.

// --- Construction and Destruction ---------------------------------------- //
MolecularGraph::MolecularGraph()
{
}

MolecularGraph::MolecularGraph(const Molecule *molecule)
    : m_molecule(molecule),
      m_atoms(molecule->atomCount()),
      m_bonds(molecule->bondCount()),
      m_adjacencyList(molecule->atomCount())
{
    for(int i = 0; i < molecule->atomCount(); i++){
        m_atoms[i] = molecule->atom(i);
    }

    for(int i = 0; i < molecule->bondCount(); i++){
        Bond *bond = molecule->bond(i);
        m_bonds[i] = bond;
        addBond(bond->atom1()->index(), bond->atom2()->index());
    }

    initializeLabels();
}

MolecularGraph::MolecularGraph(const Fragment *fragment)
    : m_molecule(fragment->molecule()),
      m_atoms(fragment->atoms()),
      m_bonds(fragment->bondCount()),
      m_adjacencyList(fragment->atomCount())
{
    const std::vector<Atom *> &atoms = fragment->atoms();

    for(int i = 0; i < fragment->bondCount(); i++){
        Bond *bond = fragment->bonds()[i];
        m_bonds[i] = bond;
        addBond(std::distance(atoms.begin(), std::find(atoms.begin(), atoms.end(), bond->atom1())),
                std::distance(atoms.begin(), std::find(atoms.begin(), atoms.end(), bond->atom2())));
    }

    initializeLabels();
}

MolecularGraph::MolecularGraph(const std::vector<Atom *> &atoms)
    : m_molecule(0),
      m_atoms(atoms),
      m_adjacencyList(atoms.size())
{
    if(atoms.empty()){
        return;
    }

    m_molecule = atoms[0]->molecule();

    for(unsigned int i = 0; i < atoms.size(); i++){
        Atom *a = atoms[i];

        for(unsigned int j = i + 1; j < atoms.size(); j++){
            Atom *b = atoms[j];
            Bond *bond = a->bondTo(b);

            if(bond){
                m_bonds.push_back(bond);
                addBond(i, j);
            }
        }
    }

    initializeLabels();
}

MolecularGraph::~MolecularGraph()
{
}

// --- Properties ---------------------------------------------------------- //
const Molecule* MolecularGraph::molecule() const
{
    return m_molecule;
}

Atom* MolecularGraph::atom(unsigned int index) const
{
    return m_atoms[index];
}

Bond* MolecularGraph::bond(unsigned int index) const
{
    return m_bonds[index];
}

unsigned int MolecularGraph::bond(unsigned int i, unsigned int j) const
{
    Atom *a = atom(i);
    Atom *b = atom(j);

    return indexOf(a->bondTo(b));
}

int MolecularGraph::indexOf(const Atom *atom) const
{
    return std::find(m_atoms.begin(), m_atoms.end(), atom) - m_atoms.begin();
}

int MolecularGraph::indexOf(const Bond *bond) const
{
    return std::find(m_bonds.begin(), m_bonds.end(), bond) - m_bonds.begin();
}

unsigned int MolecularGraph::size() const
{
    return atomCount();
}

bool MolecularGraph::isEmpty() const
{
    return m_atoms.empty();
}

unsigned int MolecularGraph::atomCount() const
{
    return m_atoms.size();
}

unsigned int MolecularGraph::bondCount() const
{
    return m_bonds.size();
}

const std::vector<int>& MolecularGraph::neighbors(unsigned int index) const
{
    return m_adjacencyList[index];
}

unsigned int MolecularGraph::neighborCount(unsigned int index) const
{
    return m_adjacencyList[index].size();
}

bool MolecularGraph::isAdjacent(unsigned int i, unsigned int j) const
{
    return std::find(neighbors(i).begin(), neighbors(i).end(), j) != neighbors(i).end();
}

// --- Labels -------------------------------------------------------------- //
void MolecularGraph::setAtomLabel(unsigned int atom, int label)
{
    m_atomLabels[atom] = label;
}

int MolecularGraph::atomLabel(unsigned int atom) const
{
    return m_atomLabels[atom];
}

void MolecularGraph::setBondLabel(unsigned int bond, int label)
{
    m_bondLabels[bond] = label;
}

int MolecularGraph::bondLabel(unsigned int bond) const
{
    return m_bondLabels[bond];
}

// --- Static Methods ------------------------------------------------------ //
MolecularGraph* MolecularGraph::cyclicGraph(const Molecule *molecule)
{
    return cyclicGraph(molecule->atoms());
}

MolecularGraph* MolecularGraph::cyclicGraph(const Fragment *fragment)
{
    return cyclicGraph(fragment->atoms());
}

MolecularGraph* MolecularGraph::cyclicGraph(const std::vector<Atom *> &atoms)
{
    std::vector<Atom *> nonterminalAtoms;

    foreach(Atom *atom, atoms){
        if(atom->neighborCount() >= 2){
            nonterminalAtoms.push_back(atom);
        }
    }

    MolecularGraph *graph = new MolecularGraph(nonterminalAtoms);
    graph->cyclicize();

    return graph;
}

MolecularGraph* MolecularGraph::hydrogenDepletedGraph(const Molecule *molecule)
{
    std::vector<Atom *> heavyAtoms;

    foreach(Atom *atom, molecule->atoms()){
        if(!atom->isTerminalHydrogen()){
            heavyAtoms.push_back(atom);
        }
    }

    MolecularGraph *graph = new MolecularGraph(heavyAtoms);
    graph->initializeLabels();
    return graph;
}

std::map<Atom *, Atom *> MolecularGraph::isomorphism(const MolecularGraph *a, const MolecularGraph *b)
{
    // if graph 'a' is larger than graph 'b' there cannot
    // be an isomorphism so just return an empty mapping
    if(a->size() > b->size()){
        return std::map<Atom *, Atom *>();
    }

    return isomorphism_vf2(a, b);
}

// --- Internal Methods ---------------------------------------------------- //
void MolecularGraph::addBond(unsigned int i, unsigned int j)
{
    m_adjacencyList[i].push_back(j);
    m_adjacencyList[j].push_back(i);
}

void MolecularGraph::removeBond(unsigned int i, unsigned int j)
{
    m_adjacencyList[i].erase(std::remove(m_adjacencyList[i].begin(), m_adjacencyList[i].end(), j));
    m_adjacencyList[j].erase(std::remove(m_adjacencyList[j].begin(), m_adjacencyList[j].end(), i));
}

void MolecularGraph::cyclicize()
{
    // remove all atoms with less than two neighbors
    bool done = false;
    unsigned int atomCount = m_atoms.size();

    while(!done){
        done = true;

        for(unsigned int i = 0; i < m_atoms.size(); i++){
            if(neighborCount(i) && neighborCount(i) < 2){
                foreach(unsigned int neighbor, neighbors(i)){
                    removeBond(i, neighbor);
                }

                atomCount--;
                done = false;
            }
        }
    }

    // compact atom list
    std::vector<Atom *> atoms;

    for(unsigned int i = 0; i < m_atoms.size(); i++){
        if(neighborCount(i) > 1){
            atoms.push_back(atom(i));
        }
    }

    // set new atoms
    m_atoms = atoms;

    // re-build adjacency list
    m_adjacencyList.clear();
    m_adjacencyList.resize(atomCount);
    m_bonds.clear();

    for(unsigned int i = 0; i < m_atoms.size(); i++){
        Atom *a = atom(i);

        for(unsigned int j = i + 1; j < m_atoms.size(); j++){
            Atom *b = atom(j);
            Bond *bond = a->bondTo(b);

            if(bond){
                m_bonds.push_back(bond);
                addBond(i, j);
            }
        }
    }
}

void MolecularGraph::initializeLabels()
{
    // set every atom label to the atom's atomic number
    m_atomLabels.resize(m_atoms.size());
    for(unsigned int i = 0; i < m_atoms.size(); i++){
        m_atomLabels[i] = m_atoms[i]->atomicNumber();
    }

    // set every bond label to the bond's order
    m_bondLabels.resize(m_bonds.size());
    for(unsigned int i = 0; i < m_bonds.size(); i++){
        m_bondLabels[i] = m_bonds[i]->order();
    }
}

std::vector<Ring *> MolecularGraph::sssr(const Molecule *molecule)
{
    std::vector<Ring *> rings;

    Q_FOREACH(const Fragment *fragment, molecule->fragments()){
        std::vector<Ring *> fragmentRings = sssr(fragment);
        rings.insert(rings.end(), fragmentRings.begin(), fragmentRings.end());
    }

    return rings;
}

std::vector<Ring *> MolecularGraph::sssr(const Fragment *fragment)
{
    MolecularGraph *graph = cyclicGraph(fragment);
    std::vector<Ring *> rings = sssr_rpPath(graph);
    delete graph;
    return rings;
}

} // end chemkit namespace
