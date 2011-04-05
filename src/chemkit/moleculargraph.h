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

#ifndef CHEMKIT_MOLECULARGRAPH_H
#define CHEMKIT_MOLECULARGRAPH_H

#include "chemkit.h"

#include <vector>

namespace chemkit {

class Atom;
class Bond;
class Ring;
class Molecule;
class Fragment;
class AtomMapping;

class CHEMKIT_EXPORT MolecularGraph
{
    public:
        // construction and destruction
        MolecularGraph(const Molecule *molecule);
        MolecularGraph(const Fragment *fragment);
        MolecularGraph(const std::vector<Atom *> &atoms);
        ~MolecularGraph();

        // properties
        const Molecule* molecule() const;
        Atom* atom(unsigned int index) const;
        Bond* bond(unsigned int index) const;
        unsigned int bond(unsigned int i, unsigned int j) const;
        int indexOf(const Atom *atom) const;
        int indexOf(const Bond *bond) const;
        unsigned int size() const;
        bool isEmpty() const;
        unsigned int atomCount() const;
        unsigned int bondCount() const;
        const std::vector<int>& neighbors(unsigned int index) const;
        unsigned int neighborCount(unsigned int index) const;
        bool isAdjacent(unsigned int i, unsigned int j) const;

        // labels
        void setAtomLabel(unsigned int atom, int label);
        int atomLabel(unsigned int atom) const;
        void setBondLabel(unsigned int bond, int label);
        int bondLabel(unsigned int bond) const;

        // static methods
        static MolecularGraph* cyclicGraph(const Molecule *molecule);
        static MolecularGraph* cyclicGraph(const Fragment *fragment);
        static MolecularGraph* cyclicGraph(const std::vector<Atom *> &atoms);
        static MolecularGraph* hydrogenDepletedGraph(const Molecule *molecule);
        static AtomMapping isomorphism(const MolecularGraph *a, const MolecularGraph *b);

    private:
        MolecularGraph();
        void addBond(unsigned int i, unsigned int j);
        void removeBond(unsigned int i, unsigned int j);
        void cyclicize();
        void initializeLabels();
        static std::vector<Ring *> sssr(const Molecule *molecule);
        static std::vector<Ring *> sssr(const Fragment *fragment);
        static std::vector<Ring *> sssr_rpPath(const MolecularGraph *graph);
        static AtomMapping isomorphism_vf2(const MolecularGraph *a, const MolecularGraph *b);

        friend class Molecule;

    private:
        const Molecule *m_molecule;
        std::vector<Atom *> m_atoms;
        std::vector<Bond *> m_bonds;
        std::vector<std::vector<int> > m_adjacencyList;
        std::vector<int> m_atomLabels;
        std::vector<int> m_bondLabels;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULARGRAPH_H
