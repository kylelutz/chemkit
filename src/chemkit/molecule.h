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

#ifndef CHEMKIT_MOLECULE_H
#define CHEMKIT_MOLECULE_H

#include "chemkit.h"

#include <map>
#include <string>
#include <vector>
#include <cstddef>

#include <boost/range/iterator_range.hpp>

#include "moiety.h"
#include "point3.h"
#include "element.h"
#include "variant.h"
#include "vector3.h"

namespace chemkit {

class Atom;
class Bond;
class Ring;
class Fragment;
class CoordinateSet;
class MoleculePrivate;
class MoleculeWatcher;
class Stereochemistry;
class DiagramCoordinates;
class InternalCoordinates;
class CartesianCoordinates;

class CHEMKIT_EXPORT Molecule
{
public:
    // enumerations
    enum ChangeType {
        AtomAdded,
        AtomRemoved,
        AtomElementChanged,
        AtomMassNumberChanged,
        AtomPartialChargeChanged,
        AtomPositionChanged,
        AtomChiralityChanged,
        BondAdded,
        BondRemoved,
        BondOrderChanged,
        NameChanged
    };

    enum CompareFlag {
        CompareAtomsOnly = 0x00,
        CompareHydrogens = 0x01,
        CompareAromaticity = 0x02
    };

    // typedefs
    typedef boost::iterator_range<std::vector<Atom *>::const_iterator> AtomRange;
    typedef boost::iterator_range<std::vector<Bond *>::const_iterator> BondRange;

    // construction and destruction
    Molecule();
    Molecule(const std::string &formula, const std::string &format);
    Molecule(const Molecule &molecule);
    virtual ~Molecule();

    // properties
    void setName(const std::string &name);
    std::string name() const;
    std::string formula() const;
    std::string formula(const std::string &format) const;
    Variant descriptor(const std::string &name) const;
    inline size_t size() const;
    inline bool isEmpty() const;
    Real mass() const;
    void setData(const std::string &name, const Variant &value);
    Variant data(const std::string &name) const;

    // structure
    Atom* addAtom(const Element &element);
    Atom* addAtomCopy(const Atom *atom);
    void removeAtom(Atom *atom);
    Atom* atom(size_t index) const;
    inline std::vector<Atom *> atoms() const;
    inline size_t atomCount() const;
    size_t atomCount(const Element &element) const;
    inline AtomRange atomRange() const;
    bool contains(const Atom *atom) const;
    bool contains(const Element &element) const;
    Bond* addBond(Atom *a, Atom *b, int order = 1);
    Bond* addBond(size_t a, size_t b, int order = 1);
    void removeBond(Bond *bond);
    void removeBond(Atom *a, Atom *b);
    void removeBond(size_t a, size_t b);
    Bond* bond(size_t index) const;
    Bond* bond(const Atom *a, const Atom *b) const;
    Bond* bond(size_t a, size_t b) const;
    std::vector<Bond *> bonds() const;
    size_t bondCount() const;
    BondRange bondRange() const;
    bool contains(const Bond *bond) const;
    void clear();

    // comparison
    bool equals(const Molecule *molecule, int flags = 0) const;
    bool contains(const Molecule *molecule, int flags = 0) const;
    bool isSubstructureOf(const Molecule *molecule, int flags = 0) const;
    std::map<Atom *, Atom *> mapping(const Molecule *molecule, int flags = 0) const;
    Moiety find(const Molecule *moiety, int flags = 0) const;

    // ring perception
    Ring* ring(size_t index) const;
    std::vector<Ring *> rings() const;
    size_t ringCount() const;

    // fragment perception
    Fragment* fragment(size_t index) const;
    std::vector<Fragment *> fragments() const;
    size_t fragmentCount() const;
    bool isFragmented() const;
    void removeFragment(Fragment *fragment);

    // coordinates
    CartesianCoordinates* coordinates() const;
    void addCoordinateSet(CoordinateSet *coordinates);
    void addCoordinateSet(CartesianCoordinates *coordinates);
    void addCoordinateSet(InternalCoordinates *coordinates);
    void addCoordinateSet(DiagramCoordinates *coordinates);
    bool removeCoordinateSet(CoordinateSet *coordinates);
    bool deleteCoordinateSet(CoordinateSet *coordinates);
    CoordinateSet* coordinateSet(size_t index) const;
    std::vector<CoordinateSet *> coordinateSets() const;
    size_t coordinateSetCount() const;

    // geometry
    Real distance(const Atom *a, const Atom *b) const;
    Real bondAngle(const Atom *a, const Atom *b, const Atom *c) const;
    Real torsionAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const;
    Real wilsonAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const;
    void setCenter(const Point3 &position);
    void setCenter(Real x, Real y, Real z);
    Point3 center() const;
    Point3 centerOfMass() const;
    void moveBy(const Vector3 &vector);
    void moveBy(Real dx, Real dy, Real dz);
    void rotate(const Vector3 &axis, Real angle);

    // operators
    Molecule& operator=(const Molecule &molecule);

private:
    // internal methods
    void setRingsPerceived(bool perceived) const;
    bool ringsPerceived() const;
    void setFragmentsPerceived(bool perceived) const;
    bool fragmentsPerceived() const;
    void perceiveFragments() const;
    Fragment* fragment(const Atom *atom) const;
    void notifyWatchers(ChangeType type);
    void notifyWatchers(const Atom *atom, ChangeType type);
    void notifyWatchers(const Bond *bond, ChangeType type);
    void addWatcher(MoleculeWatcher *watcher) const;
    void removeWatcher(MoleculeWatcher *watcher) const;
    bool isSubsetOf(const Molecule *molecule, int flags) const;
    Stereochemistry* stereochemistry();

    friend class Atom;
    friend class Bond;
    friend class MoleculeWatcher;

private:
    MoleculePrivate* const d;
    std::vector<Atom *> m_atoms;
    std::vector<Element> m_elements;
    mutable CartesianCoordinates *m_coordinates;
    Stereochemistry *m_stereochemistry;
};

} // end chemkit namespace

#include "molecule-inline.h"

#endif // CHEMKIT_MOLECULE_H
