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

#ifndef CHEMKIT_NUCLEOTIDE_H
#define CHEMKIT_NUCLEOTIDE_H

#include "chemkit.h"

#include "residue.h"

namespace chemkit {

class Atom;
class Ring;
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
    void setType(const std::string &symbol);
    NucleotideType type() const;
    char letter() const CHEMKIT_OVERRIDE;
    std::string symbol() const;
    std::string name() const;
    void setSugarType(SugarType type);
    SugarType sugarType() const;
    bool isPurine() const;
    bool isPyrimidine() const;

private:
    NucleotidePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_NUCLEOTIDE_H
