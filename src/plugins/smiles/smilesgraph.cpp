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

#include "smilesgraph.h"

#include <map>
#include <set>
#include <queue>
#include <algorithm>

#include <boost/range.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

#include "smiles.h"

// === SmilesGraphNode ===================================================== //
SmilesGraphNode::SmilesGraphNode(const chemkit::Atom *atom)
    : m_atom(atom),
      m_hydrogenCount(0),
      m_parent(0),
      m_bondOrder(0)
{
}

SmilesGraphNode::~SmilesGraphNode()
{
    foreach(SmilesGraphNode *child, m_children){
        delete child;
    }
}

void SmilesGraphNode::setParent(SmilesGraphNode *parent, int bondOrder)
{
    m_parent = parent;
    m_bondOrder = bondOrder;
    parent->m_children.push_back(this);
}

SmilesGraphNode* SmilesGraphNode::parent() const
{
    return m_parent;
}

int SmilesGraphNode::childCount() const
{
    return m_children.size();
}

std::vector<SmilesGraphNode *> SmilesGraphNode::children() const
{
    return m_children;
}

void SmilesGraphNode::setHydrogenCount(int hydrogenCount)
{
    m_hydrogenCount = hydrogenCount;
}

int SmilesGraphNode::hydrogenCount() const
{
    return m_hydrogenCount;
}

void SmilesGraphNode::addRing(int ringNumber, int bondOrder)
{
    m_rings.push_back(ringNumber);
    m_ringBondOrders.push_back(bondOrder);
}

std::string SmilesGraphNode::toString(bool kekulize) const
{
    std::stringstream string;

    if(m_bondOrder == 0){
        // do nothing
    }
    else if(!kekulize && isAromaticAtom(m_atom) && isAromaticAtom(m_parent->atom())){
        // do nothing
    }
    else if(m_bondOrder == chemkit::Bond::Double){
        string << "=";
    }
    else if(m_bondOrder == chemkit::Bond::Triple){
        string << "#";
    }
    else if(m_bondOrder == chemkit::Bond::Quadruple){
        string << "$";
    }

    if(!kekulize && isAromaticAtom(m_atom)){
        if(m_atom->is(chemkit::Atom::Nitrogen) && m_atom->neighborCount(chemkit::Atom::Hydrogen) == 1){
            string << "[nH]";
        }
        else{
            string << boost::to_lower_copy(m_atom->symbol());
        }
    }
    else if(isOrganicAtom(m_atom)){
        string << m_atom->symbol().c_str();
    }
    else{
        string << "[";

        // mass number
        if(isIsotope(m_atom)){
            string << m_atom->massNumber();
        }

        string << m_atom->symbol().c_str();

        if(m_hydrogenCount > 0){
            string << "H";

            if(m_hydrogenCount > 1){
                string << m_hydrogenCount;
            }
        }

        int charge = m_atom->formalCharge();
        if(charge > 0){
            string << "+";
        }
        else if(charge < 0){
            string << "-";
        }

        if(std::abs(charge) > 1){
            string << std::abs(charge);
        }

        string << "]";
    }

    for(size_t i = 0; i < m_rings.size(); i++){
        int ringNumber = m_rings[i];
        int bondOrder = m_ringBondOrders[i];

        if(isAromaticAtom(m_atom)){
        }
        else if(bondOrder == chemkit::Bond::Double)
            string << "=";
        else if(bondOrder == chemkit::Bond::Triple)
            string << "#";

        if(ringNumber > 9){
            string << "%";
        }

        string << ringNumber;
    }

    return string.str();
}

void SmilesGraphNode::write(std::stringstream &string, bool kekulize) const
{
    string << toString(kekulize);

    if(childCount() == 1){
        m_children[0]->write(string, kekulize);
    }
    else if(childCount() > 1){
        std::vector<SmilesGraphNode *> children = m_children;
        const SmilesGraphNode *firstChild = children[0];
        children.erase(children.begin());

        foreach(SmilesGraphNode *child, children){
            string << "(";
            child->write(string, kekulize);
            string << ")";
        }

        firstChild->write(string, kekulize);
    }
}

// === SmilesGraph ========================================================= //
SmilesGraph::SmilesGraph(const chemkit::Molecule *molecule)
{
    std::set<const chemkit::Atom *> visitedAtoms;
    std::set<const chemkit::Ring *> visitedRings;

    std::multimap<const chemkit::Atom *, int> ringClosingAtoms;
    std::set<const chemkit::Bond *> ringBonds;

    // neighbor count for each atom without implicit hydrogens
    std::vector<int> neighborCounts(molecule->size());
    for(size_t i = 0; i < molecule->size(); i++){
        const chemkit::Atom *atom = molecule->atom(i);

        int neighborCount = 0;
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(!isImplicitHydrogen(neighbor)){
                neighborCount++;
            }
        }

        neighborCounts[i] = neighborCount;
    }

    while(visitedAtoms.size() != size_t(molecule->size())){
        const chemkit::Atom *rootAtom = 0;
        foreach(const chemkit::Atom *atom, molecule->atoms()){
            if(visitedAtoms.find(atom) != visitedAtoms.end()){
                continue;
            }

            if(!isImplicitHydrogen(atom)){
                rootAtom = atom;
                break;
            }
        }

        SmilesGraphNode *rootNode = new SmilesGraphNode(rootAtom);
        m_rootNodes.push_back(rootNode);
        visitedAtoms.insert(rootAtom);

        int ringNumber = 1;

        std::queue<SmilesGraphNode *> queue;
        queue.push(rootNode);

        while(!queue.empty()){
            SmilesGraphNode *parentNode = queue.front();
            queue.pop();
            const chemkit::Atom *atom = parentNode->atom();

            int hydrogenCount = 0;

            foreach(const chemkit::Ring *ring, atom->rings()){
                if(visitedRings.find(ring) != visitedRings.end()){
                    continue;
                }
                if(neighborCounts[atom->index()] <= 1){
                    break;
                }

                const chemkit::Atom *ringClosingAtom = 0;
                foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                    if(!ring->contains(neighbor)){
                        continue;
                    }
                    else if(ringBonds.find(atom->bondTo(neighbor)) != ringBonds.end()){
                        continue;
                    }
                    else if(neighborCounts[neighbor->index()] <= 1){
                        continue;
                    }

                    ringClosingAtom = neighbor;
                }

                if(!ringClosingAtom){
                    continue;
                }

                const chemkit::Bond *bond = atom->bondTo(ringClosingAtom);
                ringClosingAtoms.insert(std::make_pair(ringClosingAtom, ringNumber));
                ringBonds.insert(bond);
                parentNode->addRing(ringNumber, bond->order());

                neighborCounts[bond->atom1()->index()]--;
                neighborCounts[bond->atom2()->index()]--;
                visitedRings.insert(ring);
                ringNumber++;
            }

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(visitedAtoms.find(neighbor) != visitedAtoms.end()){
                    continue;
                }
                else if(isImplicitHydrogen(neighbor)){
                    hydrogenCount++;
                    visitedAtoms.insert(neighbor);
                    continue;
                }

                const chemkit::Bond *bond = atom->bondTo(neighbor);
                if(ringBonds.find(bond) != ringBonds.end()){
                    continue;
                }

                SmilesGraphNode *node = new SmilesGraphNode(neighbor);
                visitedAtoms.insert(neighbor);
                node->setParent(parentNode, bond->order());

                if(ringClosingAtoms.find(neighbor) != ringClosingAtoms.end()){
                    std::pair<std::multimap<const chemkit::Atom *, int>::iterator,
                              std::multimap<const chemkit::Atom *, int>::iterator> range;
                    range = ringClosingAtoms.equal_range(neighbor);

                    std::multimap<const chemkit::Atom *, int>::iterator iter;
                    for(iter = range.first; iter != range.second; ++iter){
                        node->addRing(iter->second, 0);
                    }

                    ringClosingAtoms.erase(neighbor);
                }

                queue.push(node);
            }

            parentNode->setHydrogenCount(hydrogenCount);
        }
    }
}

SmilesGraph::~SmilesGraph()
{
    foreach(SmilesGraphNode *node, m_rootNodes){
        delete node;
    }
}

std::string SmilesGraph::toString(bool kekulize) const
{
    if(m_rootNodes.empty()){
        return std::string();
    }
    else if(m_rootNodes.size() == 1){
        std::stringstream stream;
        m_rootNodes.front()->write(stream, kekulize);
        return stream.str();
    }
    else{
        std::stringstream stream;
        foreach(const SmilesGraphNode *rootNode, m_rootNodes){
            rootNode->write(stream, kekulize);
            stream << ".";
        }

        // remove trailing '.'
        std::string string = stream.str();
        string.erase(string.end()-1);

        return string;
    }
}
