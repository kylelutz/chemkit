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

#ifndef CHEMKIT_NUCLEICACID_H
#define CHEMKIT_NUCLEICACID_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class Molecule;
class Nucleotide;
class NucleicAcidChain;
class NucleicAcidPrivate;

class CHEMKIT_EXPORT NucleicAcid
{
    public:
        // construction and destruction
        NucleicAcid();
        ~NucleicAcid();

        // properties
        int size() const;
        Molecule* molecule();
        const Molecule* molecule() const;

        // structure
        NucleicAcidChain* addChain();
        void removeChain(NucleicAcidChain *chain);
        NucleicAcidChain* chain(int index);
        const NucleicAcidChain* chain(int index) const;
        QList<NucleicAcidChain *> chains();
        QList<const NucleicAcidChain *> chains() const;
        int chainCount() const;
        int residueCount() const;

    private:
        NucleicAcidPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_NUCLEICACID_H
