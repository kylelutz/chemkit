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

#ifndef CHEMKIT_RING_H
#define CHEMKIT_RING_H

#include "chemkit.h"

#include <vector>

#include <boost/range/iterator_range.hpp>

namespace chemkit {

class Atom;
class Bond;
class Element;
class Fragment;
class Molecule;

class CHEMKIT_EXPORT Ring
{
public:
    // typedefs
    typedef boost::iterator_range<std::vector<Atom *>::const_iterator> AtomRange;

    // properties
    inline size_t size() const;
    inline Molecule* molecule() const;
    inline Fragment* fragment() const;

    // structure
    inline Atom* atom(size_t index) const;
    inline std::vector<Atom *> atoms() const;
    inline size_t atomCount() const;
    size_t atomCount(const Element &element) const;
    inline AtomRange atomRange() const;
    Bond* bond(size_t index) const;
    std::vector<Bond *> bonds() const;
    size_t bondCount() const;
    std::vector<Bond *> exocyclicBonds() const;
    size_t exocyclicBondCount() const;
    inline bool contains(const Atom *atom) const;
    inline bool contains(const Bond *bond) const;
    inline bool contains(const Element &element) const;
    size_t heteroatomCount() const;
    bool isHeterocycle() const;
    Atom* root() const;
    int position(const Atom *atom, const Atom *root = 0) const;
    bool isFusedTo(const Ring *ring) const;

    // aromaticity
    bool isAromatic() const;

private:
    Ring(std::vector<Atom *> path);
    ~Ring();

    // internal methods
    bool isValid() const;
    const Atom *nextAtom(const Atom *atom) const;
    const Atom *previousAtom(const Atom *atom) const;
    const Bond *nextBond(const Atom *atom) const;
    const Bond *previousBond(const Atom *atom) const;
    bool isPlanar() const;
    size_t piElectronCount() const;

    CHEMKIT_DISABLE_COPY(Ring)

    friend class Molecule;
    friend class MolecularGraph;

private:
    std::vector<Atom *> m_atoms;
};

} // end chemkit namespace

#include "ring-inline.h"

#endif // CHEMKIT_RING_H
