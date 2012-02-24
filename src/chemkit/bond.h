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

#ifndef CHEMKIT_BOND_H
#define CHEMKIT_BOND_H

#include "chemkit.h"

#include <vector>

#include <boost/function.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include "point3.h"
#include "vector3.h"
#include "stereochemistry.h"

namespace chemkit {

class Atom;
class Ring;
class Element;
class Fragment;
class Molecule;

class CHEMKIT_EXPORT Bond
{
public:
    // typedefs
    typedef unsigned char BondOrderType;
    typedef boost::iterator_range<
                boost::filter_iterator<
                    boost::function<bool (const Ring *)>,
                    boost::iterator_range<
                        std::vector<Ring *>::const_iterator>::const_iterator> > RingRange;

    // enumerations
    enum BondType{
        Single = 1,
        Double = 2,
        Triple = 3,
        Quadruple = 4
    };

    // properties
    Atom* atom(size_t index) const;
    Atom* atom1() const;
    Atom* atom2() const;
    Atom* otherAtom(const Atom *atom) const;
    void setOrder(BondOrderType order);
    BondOrderType order() const;
    bool is(BondOrderType order) const;
    Real polarity() const;
    Vector3 dipoleMoment() const;
    inline Molecule* molecule() const;
    Fragment* fragment() const;
    inline size_t index() const;

    // structure
    bool contains(const Atom *atom) const;
    bool contains(const Element &element) const;
    bool containsBoth(const Atom *a, const Atom *b) const;
    bool containsBoth(const Element &a, const Element &b) const;
    bool isTerminal() const;

    // ring perception
    Ring* ring(size_t index) const;
    RingRange rings() const;
    size_t ringCount() const;
    bool isInRing() const;
    bool isInRing(size_t size) const;
    Ring* smallestRing() const;
    bool isAromatic() const;

    // geometry
    Point3 center() const;
    Real length() const;

    // stereochemistry
    void setStereochemistry(Stereochemistry::Type stereochemistry);
    Stereochemistry::Type stereochemistry() const;

private:
    Bond(Molecule *molecule, size_t index);
    ~Bond();

    CHEMKIT_DISABLE_COPY(Bond)

    friend class Molecule;

private:
    Molecule *m_molecule;
    size_t m_index;
};

} // end chemkit namespace

#include "bond-inline.h"

#endif // CHEMKIT_BOND_H
