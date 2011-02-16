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

#ifndef CHEMKIT_NUCLEOTIDE_H
#define CHEMKIT_NUCLEOTIDE_H

#include "chemkit.h"

#include "residue.h"

namespace chemkit {

class Atom;
class Ring;
class NucleicAcid;
class NucleicAcidChain;
class NucleotidePrivate;

class CHEMKIT_EXPORT Nucleotide : public Residue
{
    public:
        // enumerations
        enum NucleotideType{
            Adenine,
            Guanine,
            Cytosine,
            Thymine,
            Uracil,
            UnspecifiedType
        };

        enum SugarType{
            Ribose,
            Deoxyribose
        };

        // construction and destruction
        Nucleotide(Molecule *molecule);
        ~Nucleotide();

        // properties
        void setType(NucleotideType type);
        void setType(const QString &symbol);
        NucleotideType type() const;
        QString letter() const;
        QString symbol() const;
        QString name() const;
        void setSugarType(SugarType type);
        SugarType sugarType() const;
        bool isPurine() const;
        bool isPyrimidine() const;
        NucleicAcid* nucleicAcid() const;
        NucleicAcidChain* nucleicAcidChain() const;

    private:
        NucleotidePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_NUCLEOTIDE_H
