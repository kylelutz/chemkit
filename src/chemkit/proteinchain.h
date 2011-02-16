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

#ifndef CHEMKIT_PROTEINCHAIN_H
#define CHEMKIT_PROTEINCHAIN_H

#include "chemkit.h"

#include <QtCore>

#include "aminoacid.h"

namespace chemkit {

class Protein;
class Molecule;
class AminoAcid;
class ProteinChainPrivate;

class CHEMKIT_EXPORT ProteinChain
{
    public:
        // properties
        int size() const;
        Protein* protein() const;
        Molecule* molecule() const;

        // structure
        void addResidue(AminoAcid *residue);
        AminoAcid* addNewResidue();
        void removeResidue(AminoAcid *residue);
        QList<AminoAcid *> residues() const;
        AminoAcid* residue(int index) const;
        int residueCount() const;
        int residueCount(AminoAcid::AminoAcidType type) const;
        QString sequenceString() const;
        int sequenceNumber(const AminoAcid *residue) const;

    private:
        ProteinChain(Protein *protein);
        ~ProteinChain();

        friend class Protein;

    private:
        ProteinChainPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PROTEINCHAIN_H
