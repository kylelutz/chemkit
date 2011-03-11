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

#include "nucleotide.h"

namespace chemkit {

// === NucleotidePrivate =================================================== //
class NucleotidePrivate
{
    public:
        Nucleotide::NucleotideType type;
        Nucleotide::SugarType sugarType;
};

// === Nucleotide ========================================================== //
/// \class Nucleotide nucleotide.h chemkit/nucleotide.h
/// \ingroup chemkit
/// \brief The Nucleotide class represents a single nucleotide
///        residue.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new nucleotide residue.
Nucleotide::Nucleotide(Molecule *molecule)
    : Residue(molecule, Residue::NucleotideResidue),
      d(new NucleotidePrivate)
{
    d->type = UnspecifiedType;
    d->sugarType = Deoxyribose;
}

/// Destroys the nucleotide object.
Nucleotide::~Nucleotide()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the single letter symbol of the nucleotide (e.g. "G" or
/// "A").
QString Nucleotide::letter() const
{
    switch(d->type){
        case Adenine: return "A";
        case Guanine: return "G";
        case Cytosine: return "C";
        case Thymine: return "T";
        case Uracil: return "U";

        case UnspecifiedType:
        default:
            return "X";
    }
}

/// Returns the single letter symbol of the nucleotide (e.g. "G" or
/// "A"). Same as letter().
QString Nucleotide::symbol() const
{
    return letter();
}

/// Returns the name of the nucleotide (e.g. "Guanine" or "Adenine").
QString Nucleotide::name() const
{
    switch(d->type){
        case Adenine: return "Adenine";
        case Guanine: return "Guanine";
        case Cytosine: return "Cytosine";
        case Thymine: return "Thymine";
        case Uracil: return "Uracil";

        case UnspecifiedType:
        default:
            return "Unspecified";
    }
}

/// Sets the nucleotide type.
void Nucleotide::setType(Nucleotide::NucleotideType type)
{
    d->type = type;
}

/// Sets the nucleotide type from its one letter symbol.
void Nucleotide::setType(const QString &letter)
{
    if(letter == "A")
        setType(Adenine);
    else if(letter == "G")
        setType(Guanine);
    else if(letter == "C")
        setType(Cytosine);
    else if(letter == "T")
        setType(Thymine);
    else if(letter == "U")
        setType(Uracil);
    else if(letter == "X")
        setType(UnspecifiedType);
}

/// Returns the type of the nucleotide.
Nucleotide::NucleotideType Nucleotide::type() const
{
    return d->type;
}

void Nucleotide::setSugarType(SugarType type)
{
    d->sugarType = type;
}

Nucleotide::SugarType Nucleotide::sugarType() const
{
    return d->sugarType;
}

bool Nucleotide::isPurine() const
{
    return type() == Adenine || type() == Guanine;
}

bool Nucleotide::isPyrimidine() const
{
    return type() == Cytosine || type() == Thymine || type() == Uracil;
}

} // end chemkit namespace
