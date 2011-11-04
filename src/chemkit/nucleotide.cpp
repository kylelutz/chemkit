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
/// Returns the single letter symbol of the nucleotide (e.g. 'G' or
/// 'A').
char Nucleotide::letter() const
{
    switch(d->type){
        case Adenine: return 'A';
        case Guanine: return 'G';
        case Cytosine: return 'C';
        case Thymine: return 'T';
        case Uracil: return 'U';

        case UnspecifiedType:
        default:
            return 'X';
    }
}

/// Returns the single letter symbol of the nucleotide (e.g. "G" or
/// "A"). Same as letter().
std::string Nucleotide::symbol() const
{
    return std::string(1, letter());
}

/// Returns the name of the nucleotide (e.g. "Guanine" or "Adenine").
std::string Nucleotide::name() const
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
void Nucleotide::setType(const std::string &symbol)
{
    if(symbol == "A")
        setType(Adenine);
    else if(symbol == "G")
        setType(Guanine);
    else if(symbol == "C")
        setType(Cytosine);
    else if(symbol == "T")
        setType(Thymine);
    else if(symbol == "U")
        setType(Uracil);
    else if(symbol == "X")
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
