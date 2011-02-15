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

#ifndef CHEMKIT_PROTEIN_H
#define CHEMKIT_PROTEIN_H

#include "chemkit.h"

#include "molecule.h"
#include "aminoacid.h"
#include "proteinchain.h"

namespace chemkit {

class ProteinPrivate;

class CHEMKIT_EXPORT Protein
{
    public:
        // construction and destruction
        Protein();
        ~Protein();

        // properties
        int size() const;
        Molecule* molecule();
        const Molecule* molecule() const;

        // structure
        ProteinChain* addChain();
        void removeChain(ProteinChain *chain);
        QList<ProteinChain *> chains();
        QList<const ProteinChain *> chains() const;
        ProteinChain* chain(int index);
        const ProteinChain* chain(int index) const;
        int chainCount() const;
        QList<const AminoAcid *> residues() const;
        int residueCount() const;

    private:
        ProteinPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PROTEIN_H
