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

#ifndef CHEMKIT_NUCLEICACIDCHAIN_H
#define CHEMKIT_NUCLEICACIDCHAIN_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class Nucleotide;
class NucleicAcid;
class NucleicAcidChainPrivate;

class CHEMKIT_EXPORT NucleicAcidChain
{
    public:
        // properties
        NucleicAcid* nucleicAcid();
        const NucleicAcid* nucleicAcid() const;

        // structure
        void addResidue(Nucleotide *residue);
        void removeResidue(Nucleotide *residue);
        QList<Nucleotide *> residues();
        QList<const Nucleotide *> residues() const;
        int residueCount() const;
        bool contains(const Nucleotide *residue) const;
        QString sequenceString() const;

    private:
        NucleicAcidChain(NucleicAcid *nucleicAcid);
        ~NucleicAcidChain();

        friend class NucleicAcid;

    private:
        NucleicAcidChainPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_NUCLEICACIDCHAIN_H
