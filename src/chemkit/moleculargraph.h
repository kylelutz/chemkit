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

#include <QtCore>

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
        MolecularGraph(const QList<Atom *> &atoms);
        MolecularGraph(const QList<const Atom *> &atoms);
        ~MolecularGraph();

        // properties
        const Molecule* molecule() const { return m_molecule; }
        const Atom* atom(int index) const { return m_atoms[index]; }
        const Bond* bond(int index) const { return m_bonds[index]; }
        int bond(int i, int j) const;
        int indexOf(const Atom *atom) const { return m_atoms.indexOf(atom); }
        int indexOf(const Bond *bond) const { return m_bonds.indexOf(bond); }
        void setAtomLabel(int atom, int label);
        int atomLabel(int atom) const;
        void setBondLabel(int bond, int label);
        int bondLabel(int bond) const;
        int size() const { return atomCount(); }
        bool isEmpty() const { return size() == 0; }
        int atomCount() const { return m_atoms.size(); }
        int bondCount() const { return m_bonds.size(); }
        const QVector<int>& neighbors(int index) const;
        int neighborCount(int index) const;
        bool adjacent(int i, int j) const;

        // static methods
        static MolecularGraph* cyclicGraph(const Molecule *molecule);
        static MolecularGraph* cyclicGraph(const Fragment *fragment);
        static MolecularGraph* cyclicGraph(const QList<Atom *> &atoms);
        static MolecularGraph* cyclicGraph(const QList<const Atom *> &atoms);
        static MolecularGraph* hydrogenDepletedGraph(const Molecule *molecule);
        static MolecularGraph* hydrogenDepletedGraph(const Fragment *fragment);
        static MolecularGraph* hydrogenDepletedGraph(const QList<Atom *> &atoms);
        static MolecularGraph* hydrogenDepletedGraph(const QList<const Atom *> &atoms);
        static AtomMapping isomorphism(const MolecularGraph *a, const MolecularGraph *b);

    private:
        MolecularGraph();
        void cyclicize();
        void initializeLabels();
        static QList<Ring *> sssr(const Molecule *molecule);
        static QList<Ring *> sssr(const Fragment *fragment);
        static QList<Ring *> sssr_rpPath(const MolecularGraph *graph);
        static AtomMapping isomorphism_vf2(const MolecularGraph *a, const MolecularGraph *b);

        friend class Molecule;

    private:
        const Molecule *m_molecule;
        QVector<const Atom *> m_atoms;
        QVector<const Bond *> m_bonds;
        QVector<QVector<int> > m_adjacencyList;
        QVector<int> m_atomLabels;
        QVector<int> m_bondLabels;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULARGRAPH_H
