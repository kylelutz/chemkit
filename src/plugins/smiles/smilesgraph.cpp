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

#include "smilesgraph.h"

#include <algorithm>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>

namespace {

// Returns true if the element is a member of the organic subset.
bool isOrganicElement(const chemkit::Atom *atom)
{
    switch(atom->atomicNumber()){
        case chemkit::Atom::Boron:
        case chemkit::Atom::Carbon:
        case chemkit::Atom::Nitrogen:
        case chemkit::Atom::Oxygen:
        case chemkit::Atom::Phosphorus:
        case chemkit::Atom::Sulfur:
        case chemkit::Atom::Chlorine:
        case chemkit::Atom::Bromine:
        case chemkit::Atom::Iodine:
            return true;

        default:
            return false;
    }
}

bool isOrganicAtom(const chemkit::Atom *atom)
{
    return isOrganicElement(atom) && atom->formalCharge() == 0;
}

bool isPlanarAtom(const chemkit::Atom *atom)
{
    if(atom->is(chemkit::Atom::Carbon) && atom->neighborCount() != 3){
        return false;
    }
    else if(atom->is(chemkit::Atom::Oxygen) && atom->neighborCount() != 2){
        return false;
    }
    else if(atom->is(chemkit::Atom::Sulfur) && atom->neighborCount() != 2){
        return false;
    }

    return true;
}

// Returns true if the atom is a member of the aromatic subset.
bool isAromaticElement(const chemkit::Atom *atom)
{
    switch(atom->atomicNumber()){
        case chemkit::Atom::Boron:
        case chemkit::Atom::Carbon:
        case chemkit::Atom::Nitrogen:
        case chemkit::Atom::Oxygen:
        case chemkit::Atom::Phosphorus:
        case chemkit::Atom::Sulfur:
        case chemkit::Atom::Arsenic:
        case chemkit::Atom::Selenium:
            return true;

        default:
            return false;
    }
}

bool isAromaticRing(const chemkit::Ring *ring)
{
    foreach(const chemkit::Atom *atom, ring->atoms()){
        if(!isAromaticElement(atom)){
            return false;
        }
        else if(!isPlanarAtom(atom)){
            return false;
        }

        foreach(const chemkit::Bond *bond, atom->bonds()){
            if(ring->contains(bond)){
                continue;
            }

            if(bond->order() == chemkit::Bond::Double && !bond->otherAtom(atom)->isInRing()){
                return false;
            }
        }
    }

    return true;
}

bool isAromaticAtom(const chemkit::Atom *atom)
{
    if(!isAromaticElement(atom)){
        return false;
    }

    foreach(const chemkit::Ring *ring, atom->rings()){
        if(isAromaticRing(ring)){
            return true;
        }
    }

    return false;
}

bool isIsotope(const chemkit::Atom *atom)
{
    if(atom->is(chemkit::Atom::Hydrogen))
        return atom->massNumber() != 1;
    else
        return atom->massNumber() != atom->atomicNumber() * 2;
}

// Returns true if the atom is an implicit hydrogen atom.
bool isImplicitHydrogen(const chemkit::Atom *atom)
{
    return atom->isTerminalHydrogen() &&
           atom->massNumber() == 1 &&
           !atom->isBondedTo(chemkit::Atom::Hydrogen);
}

} // end anonymous namespace

// === SmilesGraphNode ===================================================== //
SmilesGraphNode::SmilesGraphNode(const chemkit::Atom *atom)
    : m_atom(atom),
      m_hydrogenCount(0),
      m_parent(0),
      m_bondOrder(0)
{
}

void SmilesGraphNode::setParent(SmilesGraphNode *parent, int bondOrder)
{
    m_parent = parent;
    m_bondOrder = bondOrder;
    parent->m_children.append(this);
}

SmilesGraphNode* SmilesGraphNode::parent() const
{
    return m_parent;
}

int SmilesGraphNode::childCount() const
{
    return m_children.size();
}

QList<SmilesGraphNode *> SmilesGraphNode::children() const
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
    m_rings.append(ringNumber);
    m_ringBondOrders.append(bondOrder);
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
            string << QString(m_atom->symbol().c_str()).toLower().toStdString();
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

        if(qAbs(charge) > 1){
            string << qAbs(charge);
        }

        string << "]";
    }

    for(int i = 0; i < m_rings.size(); i++){
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
        QList<SmilesGraphNode *> children = m_children;
        const SmilesGraphNode *firstChild = children.takeFirst();

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
    QSet<const chemkit::Atom *> visitedAtoms;
    QSet<const chemkit::Ring *> visitedRings;

    QHash<const chemkit::Atom *, int> ringClosingAtoms;
    QSet<const chemkit::Bond *> ringBonds;

    // neighbor count for each atom without implicit hydrogens
    std::vector<int> neighborCounts(molecule->size());
    for(int i = 0; i < molecule->size(); i++){
        const chemkit::Atom *atom = molecule->atom(i);

        int neighborCount = 0;
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            if(!isImplicitHydrogen(neighbor)){
                neighborCount++;
            }
        }

        neighborCounts[i] = neighborCount;
    }

    while(visitedAtoms.size() != molecule->size()){
        const chemkit::Atom *rootAtom = 0;
        foreach(const chemkit::Atom *atom, molecule->atoms()){
            if(visitedAtoms.contains(atom)){
                continue;
            }

            if(!isImplicitHydrogen(atom)){
                rootAtom = atom;
                break;
            }
        }

        SmilesGraphNode *rootNode = new SmilesGraphNode(rootAtom);
        m_rootNodes.append(rootNode);
        visitedAtoms.insert(rootAtom);

        int ringNumber = 1;

        QQueue<SmilesGraphNode *> queue;
        queue.enqueue(rootNode);

        while(!queue.isEmpty()){
            SmilesGraphNode *parentNode = queue.dequeue();
            const chemkit::Atom *atom = parentNode->atom();

            int hydrogenCount = 0;

            foreach(const chemkit::Ring *ring, atom->rings()){
                if(visitedRings.contains(ring)){
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
                    else if(ringBonds.contains(atom->bondTo(neighbor))){
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
                ringClosingAtoms.insertMulti(ringClosingAtom, ringNumber);
                ringBonds.insert(bond);
                parentNode->addRing(ringNumber, bond->order());

                neighborCounts[bond->atom1()->index()]--;
                neighborCounts[bond->atom2()->index()]--;
                visitedRings.insert(ring);
                ringNumber++;
            }

            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                if(visitedAtoms.contains(neighbor)){
                    continue;
                }
                else if(isImplicitHydrogen(neighbor)){
                    hydrogenCount++;
                    visitedAtoms.insert(neighbor);
                    continue;
                }

                const chemkit::Bond *bond = atom->bondTo(neighbor);
                if(ringBonds.contains(bond)){
                    continue;
                }

                SmilesGraphNode *node = new SmilesGraphNode(neighbor);
                visitedAtoms.insert(neighbor);
                node->setParent(parentNode, bond->order());

                if(ringClosingAtoms.contains(neighbor)){
                    foreach(int ring, ringClosingAtoms.values(neighbor)){
                        node->addRing(ring, 0);
                    }

                    ringClosingAtoms.remove(neighbor);
                }

                queue.enqueue(node);
            }

            parentNode->setHydrogenCount(hydrogenCount);
        }
    }
}

std::string SmilesGraph::toString(bool kekulize) const
{
    if(m_rootNodes.isEmpty()){
        return std::string();
    }
    else if(m_rootNodes.size() == 1){
        std::stringstream stream;
        m_rootNodes.first()->write(stream, kekulize);
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
