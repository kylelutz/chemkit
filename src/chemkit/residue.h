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

#ifndef CHEMKIT_RESIDUE_H
#define CHEMKIT_RESIDUE_H

#include "chemkit.h"

#include <string>

#include <QtCore>

namespace chemkit {

class Atom;
class Bond;
class Molecule;
class ResiduePrivate;

class CHEMKIT_EXPORT Residue
{
    public:
        // enumerations
        enum ResidueType{
            AminoAcidResidue,
            NucleotideResidue,
            CustomResidue
        };

        // construction and destruction
        Residue(Molecule *molecule, int type = CustomResidue);
        virtual ~Residue();

        // properties
        int residueType() const;
        int size() const;
        virtual char letter() const;
        Molecule* molecule() const;

        // structure
        void addAtom(Atom *atom);
        void removeAtom(Atom *atom);
        QList<Atom *> atoms() const;
        int atomCount() const;
        QList<Bond *> bonds() const;
        int bondCount() const;
        bool contains(const Atom *atom) const;
        bool contains(const Bond *bond) const;

        // atom types
        void setAtomType(const Atom *atom, const std::string &type);
        std::string atomType(const Atom *atom) const;
        Atom* atom(const std::string &type) const;

    private:
        Q_DISABLE_COPY(Residue)

        friend class Molecule;

    private:
        ResiduePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_RESIDUE_H
