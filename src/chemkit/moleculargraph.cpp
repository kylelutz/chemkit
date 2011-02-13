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
        const Atom *atom = molecule->atom(i);
        m_atoms[i] = atom;
    }

    for(int i = 0; i < molecule->bondCount(); i++){
        const Bond *bond = molecule->bond(i);
        m_bonds[i] = bond;
        m_adjacencyList[bond->atom1()->index()].append(bond->atom2()->index());
        m_adjacencyList[bond->atom2()->index()].append(bond->atom1()->index());
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
        const Atom *atom = fragment->atoms()[i];
        m_atoms[i] = atom;
    }

    for(int i = 0; i < fragment->bondCount(); i++){
        const Bond *bond = fragment->bonds()[i];
        m_bonds[i] = bond;
        m_adjacencyList[bond->atom1()->index()].append(bond->atom2()->index());
        m_adjacencyList[bond->atom2()->index()].append(bond->atom1()->index());
    }

    initializeLabels();
}

MolecularGraph::MolecularGraph(const QList<Atom *> &atoms)
    : m_molecule(0),
      m_atoms(atoms.size()),
      m_adjacencyList(atoms.size())
{
    if(atoms.isEmpty()){
        return;
    }

    m_molecule = atoms[0]->molecule();

    for(int i = 0; i < atoms.size(); i++){
        const Atom *atom = atoms[i];
        m_atoms[i] = atom;
    }

    for(int i = 0; i < m_atoms.size(); i++){
        const Atom *a = m_atoms[i];

        for(int j = i + 1; j < m_atoms.size(); j++){
            const Atom *b = m_atoms[j];

            const Bond *bond = a->bondTo(b);
            if(bond){
                m_bonds.append(bond);
                m_adjacencyList[i].append(j);
                m_adjacencyList[j].append(i);
            }
        }
    }

    initializeLabels();
}

MolecularGraph::MolecularGraph(const QList<const Atom *> &atoms)
    : m_molecule(0),
      m_atoms(atoms.size()),
      m_adjacencyList(atoms.size())
{
    if(atoms.isEmpty()){
        return;
    }

    m_molecule = atoms[0]->molecule();

    for(int i = 0; i < atoms.size(); i++){
        const Atom *atom = atoms[i];
        m_atoms[i] = atom;
    }

    for(int i = 0; i < m_atoms.size(); i++){
        const Atom *a = m_atoms[i];

        for(int j = i + 1; j < m_atoms.size(); j++){
            const Atom *b = m_atoms[j];

            const Bond *bond = a->bondTo(b);
            if(bond){
                m_bonds.append(bond);
                m_adjacencyList[i].append(j);
                m_adjacencyList[j].append(i);
            }
        }
    }

    initializeLabels();
}

MolecularGraph::~MolecularGraph()
{
}

// --- Properties ---------------------------------------------------------- //
int MolecularGraph::bond(int i, int j) const
{
    const Atom *a = atom(i);
    const Atom *b = atom(j);

    return indexOf(a->bondTo(b));
}

void MolecularGraph::setAtomLabel(int atom, int label)
{
    m_atomLabels[atom] = label;
}

int MolecularGraph::atomLabel(int atom) const
{
    return m_atomLabels[atom];
}

void MolecularGraph::setBondLabel(int bond, int label)
{
    m_bondLabels[bond] = label;
}

int MolecularGraph::bondLabel(int bond) const
{
    return m_bondLabels[bond];
}

int MolecularGraph::neighborCount(int index) const
{
    return m_adjacencyList[index].size();
}

const QVector<int>& MolecularGraph::neighbors(int index) const
{
    return m_adjacencyList[index];
}

bool MolecularGraph::adjacent(int i, int j) const
{
    return neighbors(i).contains(j);
}

// --- Ring Perception ----------------------------------------------------- //
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

// --- Internal Methods ---------------------------------------------------- //
void MolecularGraph::cyclicize()
{
    // remove all atoms with less than two neighbors
    int atomCount = m_atoms.size();
    bool done = false;

    while(!done){
        done = true;

        for(int i = 0; i < m_atoms.size(); i++){
            if(neighborCount(i) && neighborCount(i) < 2){
                foreach(int neighbor, neighbors(i)){
                    m_adjacencyList[neighbor].remove(m_adjacencyList[neighbor].indexOf(i));
                }

                m_adjacencyList[i].clear();
                atomCount--;

                done = false;
            }
        }
    }

    // compact atom list
    QVector<const Atom *> atoms;

    for(int i = 0; i < m_atoms.size(); i++){
        if(neighborCount(i) > 1){
            atoms.append(atom(i));
        }
    }

    m_atoms = atoms;

    // re-build adjacency list
    QVector<const Bond *> bonds;
    QVector<QVector<int> > adjacencyList(atomCount);

    for(int i = 0; i < m_atoms.size(); i++){
        const Atom *a = atom(i);

        for(int j = i + 1; j < m_atoms.size(); j++){
            const Atom *b = atom(j);
            const Bond *bond = a->bondTo(b);

            if(bond){
                bonds.append(bond);
                adjacencyList[i].append(j);
                adjacencyList[j].append(i);
            }
        }
    }

    m_bonds = bonds;
    m_adjacencyList = adjacencyList;
}

void MolecularGraph::initializeLabels()
{
    // set every atom label to the atom's atomic number
    m_atomLabels.resize(m_atoms.size());
    for(int i = 0; i < m_atoms.size(); i++){
        m_atomLabels[i] = m_atoms[i]->atomicNumber();
    }

    // set every bond label to the bond's order
    m_bondLabels.resize(m_bonds.size());
    for(int i = 0; i < m_bonds.size(); i++){
        m_bondLabels[i] = m_bonds[i]->order();
    }
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
    QList<const Atom *> constAtoms;
    foreach(const Atom *atom, atoms){
        constAtoms.append(atom);
    }

    return cyclicGraph(constAtoms);
}

MolecularGraph* MolecularGraph::cyclicGraph(const QList<const Atom *> &atoms)
{
    QList<const Atom *> nonTerminalAtoms;
    foreach(const Atom *atom, atoms){
        if(atom->neighborCount() >= 2){
            nonTerminalAtoms.append(atom);
        }
    }

    MolecularGraph *graph = new MolecularGraph(nonTerminalAtoms);
    graph->cyclicize();
    graph->initializeLabels();
    return graph;
}

MolecularGraph* MolecularGraph::hydrogenDepletedGraph(const Molecule *molecule)
{
    return hydrogenDepletedGraph(molecule->atoms());
}

MolecularGraph* MolecularGraph::hydrogenDepletedGraph(const Fragment *fragment)
{
    return hydrogenDepletedGraph(fragment->atoms());
}

MolecularGraph* MolecularGraph::hydrogenDepletedGraph(const QList<Atom *> &atoms)
{
    QList<const Atom *> constAtoms;
    foreach(const Atom *atom, atoms){
        constAtoms.append(atom);
    }

    return hydrogenDepletedGraph(constAtoms);
}

MolecularGraph* MolecularGraph::hydrogenDepletedGraph(const QList<const Atom *> &atoms)
{
    QList<const Atom *> heavyAtoms;
    foreach(const Atom *atom, atoms){
        if(!atom->isTerminalHydrogen()){
            heavyAtoms.append(atom);
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

} // end chemkit namespace
