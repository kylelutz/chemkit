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

#ifndef SMILESGRAPH_H
#define SMILESGRAPH_H

#include <QtCore>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>

class SmilesGraphNode
{
    public:
        SmilesGraphNode(const chemkit::Atom *atom);

        const chemkit::Atom* atom() const { return m_atom; }

        void setParent(SmilesGraphNode *parent, int bondOrder);
        SmilesGraphNode* parent() const;
        int childCount() const;
        QList<SmilesGraphNode *> children() const;
        QString toString(bool kekulize) const;
        void write(QString &string, bool kekulize) const;

        void setHydrogenCount(int hydrogenCount);
        int hydrogenCount() const;

        void addRing(int ringNumber, int bondOrder);

    private:
        const chemkit::Atom *m_atom;
        int m_hydrogenCount;
        SmilesGraphNode *m_parent;
        int m_bondOrder;
        QList<SmilesGraphNode *> m_children;
        QList<int> m_rings;
        QList<int> m_ringBondOrders;
};

class SmilesGraph
{
    public:
        SmilesGraph(const chemkit::Molecule *molecule);

        QString toString(bool kekulize) const;

    private:
        QList<SmilesGraphNode *> m_rootNodes;
};

#endif // SMILESGRAPH_H
