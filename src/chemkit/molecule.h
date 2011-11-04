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

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "moiety.h"
#include "point3.h"
#include "variant.h"
#include "vector3.h"
#include "conformer.h"

namespace chemkit {

class Fragment;
class Coordinates;
class MoleculePrivate;
class Stereochemistry;
class MoleculeObserver;
class InternalCoordinates;

class CHEMKIT_EXPORT Molecule
{
public:
    // enumerations
    enum ChangeType {
        AtomAdded,
        AtomRemoved,
        AtomAtomicNumberChanged,
        AtomMassNumberChanged,
        AtomPartialChargeChanged,
        AtomPositionChanged,
        AtomChiralityChanged,
        BondAdded,
        BondRemoved,
        BondOrderChanged,
        ConformerAdded,
        ConformerRemoved,
        ConformerChanged,
        NameChanged
    };

    enum CompareFlag {
        CompareAtomsOnly = 0x00,
        CompareHydrogens = 0x01,
        CompareAromaticity = 0x02
    };

    // construction and destruction
    Molecule();
    Molecule(const std::string &formula, const std::string &format);
    Molecule(const Molecule &molecule);
    ~Molecule();

    // properties
    void setName(const std::string &name);
    std::string name() const;
    std::string formula() const;
    std::string formula(const std::string &format) const;
    Variant descriptor(const std::string &name) const;
    inline int size() const;
    inline bool isEmpty() const;
    Real mass() const;
    void setData(const std::string &name, const Variant &value);
    Variant data(const std::string &name) const;

    // structure
    Atom* addAtom(const Element &element);
    Atom* addAtomCopy(const Atom *atom);
    void removeAtom(Atom *atom);
    Atom* atom(int index) const;
    inline std::vector<Atom *> atoms() const;
    inline int atomCount() const;
    int atomCount(const Element &element) const;
    bool contains(const Atom *atom) const;
    bool contains(const Element &element) const;
    Bond* addBond(Atom *a, Atom *b, int order = Bond::Single);
    Bond* addBond(int a, int b, int order = Bond::Single);
    void removeBond(Bond *bond);
    void removeBond(Atom *a, Atom *b);
    void removeBond(int a, int b);
    Bond* bond(int index) const;
    Bond* bond(const Atom *a, const Atom *b) const;
    Bond* bond(int a, int b) const;
    std::vector<Bond *> bonds() const;
    int bondCount() const;
    bool contains(const Bond *bond) const;
    void clear();

    // comparison
    bool equals(const Molecule *molecule, int flags = 0) const;
    bool contains(const Molecule *molecule, int flags = 0) const;
    bool isSubstructureOf(const Molecule *molecule, int flags = 0) const;
    std::map<Atom *, Atom *> mapping(const Molecule *molecule, int flags = 0) const;
    Moiety find(const Molecule *moiety, int flags = 0) const;

    // ring perception
    Ring* ring(int index) const;
    std::vector<Ring *> rings() const;
    int ringCount() const;

    // fragment perception
    Fragment* fragment(int index) const;
    std::vector<Fragment *> fragments() const;
    int fragmentCount() const;
    bool isFragmented() const;
    void removeFragment(Fragment *fragment);

    // geometry
    void setCoordinates(const Coordinates *coordinates);
    void setCoordinates(const InternalCoordinates *coordinates);
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
    bool hasCoordinates() const;
    void clearCoordinates();

    // conformers
    Conformer* addConformer();
    void removeConformer(Conformer *conformer);
    void setConformer(Conformer *conformer);
    Conformer* conformer() const;
    Conformer* conformer(int index) const;
    std::vector<Conformer *> conformers() const;
    int conformerCount() const;

    // operators
    Molecule& operator=(const Molecule &molecule);

private:
    // internal methods
    std::vector<Atom *> atomPathBetween(const Atom *a, const Atom *b) const;
    int atomCountBetween(const Atom *a, const Atom *b) const;
    int atomCountBetween(const Atom *a, const Atom *b, int maxCount) const;
    std::vector<Bond *> bondPathBetween(const Atom *a, const Atom *b) const;
    int bondCountBetween(const Atom *a, const Atom *b) const;
    int bondCountBetween(const Atom *a, const Atom *b, int maxCount) const;
    void setRingsPerceived(bool perceived) const;
    bool ringsPerceived() const;
    void setFragmentsPerceived(bool perceived) const;
    bool fragmentsPerceived() const;
    void perceiveFragments() const;
    Fragment* fragment(const Atom *atom) const;
    void notifyObservers(ChangeType type);
    void notifyObservers(const Atom *atom, ChangeType type);
    void notifyObservers(const Bond *bond, ChangeType type);
    void notifyObservers(const Conformer *conformer, ChangeType type);
    void addObserver(MoleculeObserver *observer) const;
    void removeObserver(MoleculeObserver *observer) const;
    bool isSubsetOf(const Molecule *molecule, int flags) const;
    Stereochemistry* stereochemistry();

    friend class Atom;
    friend class Bond;
    friend class MoleculeObserver;

private:
    MoleculePrivate* const d;
    std::vector<Atom *> m_atoms;
    std::vector<Element> m_elements;
    Stereochemistry *m_stereochemistry;
};

} // end chemkit namespace

#include "molecule-inline.h"

#endif // CHEMKIT_MOLECULE_H
