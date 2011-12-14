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

#include "aminoacid.h"

#include <boost/algorithm/string.hpp>

#include "atom.h"
#include "geometry.h"

namespace chemkit {

// === AminoAcidPrivate ==================================================== //
class AminoAcidPrivate
{
public:
    AminoAcid::AminoAcidType type;
    AminoAcid::Conformation conformation;
    Atom *alphaCarbon;
    Atom *carbonylCarbon;
    Atom *carbonylOxygen;
    Atom *aminoNitrogen;
};

// === AminoAcid =========================================================== //
/// \class AminoAcid aminoacid.h chemkit/aminoacid.h
/// \ingroup chemkit
/// \brief The AminoAcid class represents a single amino acid in a
///        protein.
///
/// \see Protein

/// \enum AminoAcid::Conformation
/// Provides names for the different conformations.
///   - \c Coil
///   - \c AlphaHelix
///   - \c BetaSheet

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new amino acid residue.
AminoAcid::AminoAcid(Molecule *molecule)
    : Residue(molecule, Residue::AminoAcidResidue),
      d(new AminoAcidPrivate)
{
    d->type = UnspecifiedType;
    d->conformation = Coil;
    d->alphaCarbon = 0;
    d->carbonylCarbon = 0;
    d->carbonylOxygen = 0;
    d->aminoNitrogen = 0;
}

/// Destroys the amino acid object.
AminoAcid::~AminoAcid()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the amino acid's type.
void AminoAcid::setType(AminoAcidType type)
{
    d->type = type;
}

/// Sets the amino acid's type from its single letter symbol (e.g.
/// "W", "E") or its three letter symbol (e.g. "Trp", "Glu").
void AminoAcid::setType(const std::string &letterOrSymbol)
{
    // set by 1-character letter
    if(letterOrSymbol.length() == 1){
        char letter = toupper(letterOrSymbol[0]);

        if(letter == 'A')
            setType(Alanine);
        else if(letter == 'R')
            setType(Arganine);
        else if(letter == 'N')
            setType(Asparagine);
        else if(letter == 'D')
            setType(AsparticAcid);
        else if(letter == 'C')
            setType(Cystenine);
        else if(letter == 'E')
            setType(GlutamicAcid);
        else if(letter == 'Q')
            setType(Glutamine);
        else if(letter == 'G')
            setType(Glycine);
        else if(letter == 'H')
            setType(Histadine);
        else if(letter == 'I')
            setType(Isoleucine);
        else if(letter == 'L')
            setType(Leucine);
        else if(letter == 'K')
            setType(Lysine);
        else if(letter == 'M')
            setType(Methionine);
        else if(letter == 'F')
            setType(Phenylalanine);
        else if(letter == 'P')
            setType(Proline);
        else if(letter == 'S')
            setType(Serine);
        else if(letter == 'T')
            setType(Threonine);
        else if(letter == 'W')
            setType(Tryptophan);
        else if(letter == 'Y')
            setType(Tyrosine);
        else if(letter == 'V')
            setType(Valine);
        else
            setType(UnspecifiedType);
    }

    // set by 3-character symbol
    else if(letterOrSymbol.length() == 3){
        std::string symbol = boost::algorithm::to_upper_copy(letterOrSymbol);

        if(symbol == "ALA")
            setType(Alanine);
        else if(symbol == "ARG")
            setType(Arganine);
        else if(symbol == "ASN")
            setType(Asparagine);
        else if(symbol == "ASP")
            setType(AsparticAcid);
        else if(symbol == "CYS")
            setType(Cystenine);
        else if(symbol == "GLU")
            setType(GlutamicAcid);
        else if(symbol == "GLN")
            setType(Glutamine);
        else if(symbol == "GLY")
            setType(Glycine);
        else if(symbol == "HIS")
            setType(Histadine);
        else if(symbol == "ILE")
            setType(Isoleucine);
        else if(symbol == "LEU")
            setType(Leucine);
        else if(symbol == "LYS")
            setType(Lysine);
        else if(symbol == "MET")
            setType(Methionine);
        else if(symbol == "PHE")
            setType(Phenylalanine);
        else if(symbol == "PRO")
            setType(Proline);
        else if(symbol == "SER")
            setType(Serine);
        else if(symbol == "THR")
            setType(Threonine);
        else if(symbol == "TRP")
            setType(Tryptophan);
        else if(symbol == "TYR")
            setType(Tyrosine);
        else if(symbol == "VAL")
            setType(Valine);
        else
            setType(UnspecifiedType);
    }
}

/// Returns the type of the amino acid.
AminoAcid::AminoAcidType AminoAcid::type() const
{
    return d->type;
}

/// Returns the name of the amino acid. (e.g. "Tryptophan" or
/// "Glutamic Acid").
std::string AminoAcid::name() const
{
    switch(type()){
        case Alanine: return "Alanine";
        case Arganine: return "Arganine";
        case Asparagine: return "Asparagine";
        case AsparticAcid: return "Aspartic Acid";
        case Cystenine: return "Cystenine";
        case GlutamicAcid: return "Glutamic Acid";
        case Glutamine: return "Glutamine";
        case Glycine: return "Glycine";
        case Histadine: return "Histadine";
        case Isoleucine: return "Isoleucine";
        case Leucine: return "Leucine";
        case Lysine: return "Lysine";
        case Methionine: return "Methionine";
        case Phenylalanine: return "Phenylalanine";
        case Proline: return "Proline";
        case Serine: return "Serine";
        case Threonine: return "Threonine";
        case Tryptophan: return "Tryptophan";
        case Tyrosine: return "Tyrosine";
        case Valine: return "Valine";

        case UnspecifiedType:
        default:
            return "Unspecified";
    }
}

/// Returns the three letter symbol of the amino acid. (e.g. "Trp" or
/// "Glu").
std::string AminoAcid::symbol() const
{
    switch(type()){
        case Alanine: return "Ala";
        case Arganine: return "Arg";
        case Asparagine: return "Asn";
        case AsparticAcid: return "Asp";
        case Cystenine: return "Cys";
        case GlutamicAcid: return "Glu";
        case Glutamine: return "Gln";
        case Glycine: return "Gly";
        case Histadine: return "His";
        case Isoleucine: return "Ile";
        case Leucine: return "Leu";
        case Lysine: return "Lys";
        case Methionine: return "Met";
        case Phenylalanine: return "Phe";
        case Proline: return "Pro";
        case Serine: return "Ser";
        case Threonine: return "Thr";
        case Tryptophan: return "Trp";
        case Tyrosine: return "Tyr";
        case Valine: return "Val";

        case UnspecifiedType:
        default:
            return "Xaa";
    }
}

/// Returns the single letter symbol of the amino acid. (e.g. 'W' or
/// 'E').
char AminoAcid::letter() const
{
    switch(type()){
        case Alanine: return 'A';
        case Arganine: return 'R';
        case Asparagine: return 'N';
        case AsparticAcid: return 'D';
        case Cystenine: return 'C';
        case GlutamicAcid: return 'E';
        case Glutamine: return 'Q';
        case Glycine: return 'G';
        case Histadine: return 'H';
        case Isoleucine: return 'I';
        case Leucine: return 'L';
        case Lysine: return 'K';
        case Methionine: return 'M';
        case Phenylalanine: return 'F';
        case Proline: return 'P';
        case Serine: return 'S';
        case Threonine: return 'T';
        case Tryptophan: return 'W';
        case Tyrosine: return 'Y';
        case Valine: return 'V';

        case UnspecifiedType:
        default:
            return 'X';
    }
}

/// Sets the conformation of the amino acid.
void AminoAcid::setConformation(Conformation conformation)
{
    d->conformation = conformation;
}

/// Returns the conformation of the amino acid.
AminoAcid::Conformation AminoAcid::conformation() const
{
    return d->conformation;
}

// --- Structure ----------------------------------------------------------- //
/// Sets the alpha carbon.
void AminoAcid::setAlphaCarbon(Atom *atom)
{
    d->alphaCarbon = atom;
}

/// Returns the alpha carbon.
Atom* AminoAcid::alphaCarbon() const
{
    return d->alphaCarbon;
}

/// Returns the carbonyl carbon.
Atom* AminoAcid::carbonylCarbon() const
{
    return d->carbonylCarbon;
}

/// Sets the carbonyl carbon.
void AminoAcid::setCarbonylCarbon(Atom *atom)
{
    d->carbonylCarbon = atom;
}

/// Sets the carbonyl oxygen.
void AminoAcid::setCarbonylOxygen(Atom *atom)
{
    d->carbonylOxygen = atom;
}

/// Returns the carbonyl oxygen.
Atom* AminoAcid::carbonylOxygen() const
{
    return d->carbonylOxygen;
}

/// Sets the amino nitrogen.
void AminoAcid::setAminoNitrogen(Atom *atom)
{
    d->aminoNitrogen = atom;
}

/// Returns the amino nitrogen.
Atom* AminoAcid::aminoNitrogen() const
{
    return d->aminoNitrogen;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the normal vector to the peptide plane.
Vector3 AminoAcid::peptidePlaneNormal() const
{
    if(!d->alphaCarbon || !d->carbonylCarbon || !d->carbonylOxygen)
        return Vector3::UnitY();

    return chemkit::geometry::planeNormal(d->alphaCarbon->position(),
                                          d->carbonylCarbon->position(),
                                          d->carbonylOxygen->position());
}

} // end chemkit namespace
