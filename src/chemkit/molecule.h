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

#ifndef CHEMKIT_MOLECULE_H
#define CHEMKIT_MOLECULE_H

#include "chemkit.h"

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "point.h"
#include "moiety.h"
#include "residue.h"
#include "fragment.h"
#include "conformer.h"
#include "atommapping.h"

namespace chemkit {

class Vector;
class MoleculePrivate;
class MoleculeWatcher;

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
            AtomResidueChanged,
            BondAdded,
            BondRemoved,
            BondOrderChanged,
            ResidueAdded,
            ResidueRemoved,
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
        Q_DECLARE_FLAGS(CompareFlags, CompareFlag)

        // construction and destruction
        Molecule();
        Molecule(const QString &formula, const QString &format);
        Molecule(const Molecule &molecule);
        ~Molecule();

        // properties
        void setName(const QString &name);
        QString name() const;
        bool hasName() const;
        QString formula() const;
        QString formula(const QString &format) const;
        QVariant descriptor(const QString &name) const;
        int size() const;
        bool isEmpty() const;
        Float mass() const;

        // structure
        Atom* addAtom(const Element &element);
        Atom* addAtomCopy(const Atom *atom);
        void removeAtom(Atom *atom);
        Atom* atom(int index) const;
        QList<Atom *> atoms() const;
        int atomCount() const;
        int atomCount(const Element &element) const;
        int indexOf(const Atom *atom) const;
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
        QList<Bond *> bonds() const;
        int bondCount() const;
        int indexOf(const Bond *bond) const;
        bool contains(const Bond *bond) const;
        void addResidue(Residue *residue);
        void removeResidue(Residue *residue);
        QList<Residue *> residues() const;
        int residueCount() const;
        void clear();

        // comparison
        bool equals(const Molecule *molecule, CompareFlags flags = CompareFlags()) const;
        bool contains(const Molecule *molecule, CompareFlags flags = CompareFlags()) const;
        bool isSubstructureOf(const Molecule *molecule, CompareFlags flags = CompareFlags()) const;
        AtomMapping mapping(const Molecule *molecule, CompareFlags flags = CompareFlags()) const;
        Moiety find(const Molecule *moiety, CompareFlags flags = CompareFlags()) const;

        // ring perception
        Ring* ring(int index) const;
        QList<Ring *> rings() const;
        int ringCount() const;

        // fragment perception
        Fragment* fragment(int index) const;
        QList<Fragment *> fragments() const;
        int fragmentCount() const;
        bool isFragmented() const;
        void removeFragment(Fragment *fragment);

        // geometry
        Float distance(const Atom *a, const Atom *b) const;
        Float bondAngle(const Atom *a, const Atom *b, const Atom *c) const;
        Float torsionAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const;
        Float wilsonAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const;
        void setCenter(const Point &position);
        void setCenter(Float x, Float y, Float z);
        Point center() const;
        Point centerOfMass() const;
        void moveBy(const Vector &vector);
        void moveBy(Float dx, Float dy, Float dz);
        void rotate(const Vector &axis, Float angle);
        bool hasCoordinates() const;
        void clearCoordinates();

        // conformers
        Conformer* addConformer();
        void removeConformer(Conformer *conformer);
        void setConformer(Conformer *conformer);
        Conformer* conformer() const;
        Conformer* conformer(int index) const;
        QList<Conformer *> conformers() const;
        int conformerCount() const;

        // operators
        Molecule& operator=(const Molecule &molecule);

    private:
        // internal methods
        QList<Atom *> atomPathBetween(const Atom *a, const Atom *b) const;
        int atomCountBetween(const Atom *a, const Atom *b) const;
        int atomCountBetween(const Atom *a, const Atom *b, int maxCount) const;
        QList<Bond *> bondPathBetween(const Atom *a, const Atom *b) const;
        int bondCountBetween(const Atom *a, const Atom *b) const;
        int bondCountBetween(const Atom *a, const Atom *b, int maxCount) const;
        void setRingsPerceived(bool perceived) const;
        bool ringsPerceived() const;
        void setFragmentsPerceived(bool perceived) const;
        bool fragmentsPerceived() const;
        Fragment* fragment(const Atom *atom) const;
        void notifyObservers(ChangeType type);
        void notifyObservers(const Atom *atom, ChangeType type);
        void notifyObservers(const Bond *bond, ChangeType type);
        void notifyObservers(const Residue *residue, ChangeType type);
        void notifyObservers(const Conformer *conformer, ChangeType type);
        void addWatcher(MoleculeWatcher *watcher) const;
        void removeWatcher(MoleculeWatcher *watcher) const;
        bool isSubsetOf(const Molecule *molecule, CompareFlags flags = CompareFlags()) const;

        friend class Atom;
        friend class Bond;
        friend class Residue;
        friend class MoleculeWatcher;

    private:
        MoleculePrivate* const d;
        QList<Atom *> m_atoms;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(Molecule::CompareFlags)

} // end chemkit namespace

#include "molecule-inline.h"

#endif // CHEMKIT_MOLECULE_H
