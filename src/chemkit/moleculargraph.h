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

#ifndef CHEMKIT_MOLECULARGRAPH_H
#define CHEMKIT_MOLECULARGRAPH_H

#include "chemkit.h"

#include <map>
#include <vector>

namespace chemkit {

class Atom;
class Bond;
class Ring;
class Molecule;
class Fragment;

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
        static std::map<Atom *, Atom *> isomorphism(const MolecularGraph *a, const MolecularGraph *b);

    private:
        MolecularGraph();
        void addBond(unsigned int i, unsigned int j);
        void removeBond(unsigned int i, unsigned int j);
        void cyclicize();
        void initializeLabels();
        static std::vector<Ring *> sssr(const Molecule *molecule);
        static std::vector<Ring *> sssr(const Fragment *fragment);
        static std::vector<Ring *> sssr_rpPath(const MolecularGraph *graph);
        static std::map<Atom *, Atom *> isomorphism_vf2(const MolecularGraph *a, const MolecularGraph *b);

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
