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

#include "moleculargraph.h"

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "molecule.h"
#include "atommapping.h"

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
      m_atoms(fragment->atomCount()),
      m_bonds(fragment->bondCount()),
      m_adjacencyList(fragment->atomCount())
{
    for(int i = 0; i < fragment->atomCount(); i++){
        m_atoms[i] = fragment->atom(i);
    }

    for(int i = 0; i < fragment->bondCount(); i++){
        Bond *bond = fragment->bonds()[i];
        m_bonds[i] = bond;
        addBond(fragment->atoms().indexOf(bond->atom1()), fragment->atoms().indexOf(bond->atom2()));
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

MolecularGraph* MolecularGraph::cyclicGraph(const QList<Atom *> &atoms)
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

AtomMapping MolecularGraph::isomorphism(const MolecularGraph *a, const MolecularGraph *b)
{
    // if graph 'a' is larger than graph 'b' there cannot
    // be an isomorphism so just return an empty mapping
    if(a->size() > b->size())
        return AtomMapping(a->molecule(), b->molecule());

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

QList<Ring *> MolecularGraph::sssr(const Molecule *molecule)
{
    QList<Ring *> rings;

    foreach(const Fragment *fragment, molecule->fragments()){
        rings += sssr(fragment);
    }

    return rings;
}

QList<Ring *> MolecularGraph::sssr(const Fragment *fragment)
{
    MolecularGraph *graph = cyclicGraph(fragment);
    QList<Ring *> rings = sssr_rpPath(graph);
    delete graph;
    return rings;
}

} // end chemkit namespace
