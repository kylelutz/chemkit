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

#ifndef SMILESGRAPH_H
#define SMILESGRAPH_H

#include <QtCore>

#include <string>
#include <sstream>

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
        std::string toString(bool kekulize) const;
        void write(std::stringstream &string, bool kekulize) const;

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

        std::string toString(bool kekulize) const;

    private:
        QList<SmilesGraphNode *> m_rootNodes;
};

#endif // SMILESGRAPH_H
