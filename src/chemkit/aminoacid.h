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

#ifndef CHEMKIT_AMINOACID_H
#define CHEMKIT_AMINOACID_H

#include "chemkit.h"

#include "residue.h"
#include "vector3.h"

namespace chemkit {

class AminoAcidPrivate;

class CHEMKIT_EXPORT AminoAcid : public Residue
{
public:
    enum Conformation {
        Coil,
        AlphaHelix,
        BetaSheet
    };

    enum AminoAcidType {
        Alanine,
        Arganine,
        Asparagine,
        AsparticAcid,
        Aspartate = AsparticAcid,
        Cystenine,
        Glutamine,
        GlutamicAcid,
        Glutamate = GlutamicAcid,
        Glycine,
        Histadine,
        Isoleucine,
        Leucine,
        Lysine,
        Methionine,
        Phenylalanine,
        Proline,
        Serine,
        Threonine,
        Tryptophan,
        Tyrosine,
        Valine,
        UnspecifiedType
    };

    // construction and destruction
    AminoAcid(Molecule *molecule);
    ~AminoAcid();

    // properties
    void setType(AminoAcidType type);
    void setType(const std::string &letterOrSymbol);
    AminoAcidType type() const;
    std::string name() const;
    std::string symbol() const;
    char letter() const CHEMKIT_OVERRIDE;
    void setConformation(Conformation conformation);
    Conformation conformation() const;

    // structure
    void setAlphaCarbon(Atom *atom);
    Atom* alphaCarbon() const;
    void setCarbonylCarbon(Atom *atom);
    Atom* carbonylCarbon() const;
    void setCarbonylOxygen(Atom *atom);
    Atom* carbonylOxygen() const;
    void setAminoNitrogen(Atom *atom);
    Atom* aminoNitrogen() const;

    // geometry
    Vector3 peptidePlaneNormal() const;

private:
    AminoAcidPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_AMINOACID_H
