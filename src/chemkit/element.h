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

#ifndef CHEMKIT_ELEMENT_H
#define CHEMKIT_ELEMENT_H

#include "chemkit.h"

#include <string>

namespace chemkit {

class CHEMKIT_EXPORT Element
{
public:
    // typedefs
    typedef unsigned char AtomicNumberType;

    // construction and destruction
    inline Element(AtomicNumberType atomicNumber = 0);
    Element(const char *symbol);
    Element(const std::string &symbol);

    // properties
    inline void setAtomicNumber(AtomicNumberType atomicNumber);
    inline AtomicNumberType atomicNumber() const;
    std::string symbol() const;
    std::string name() const;
    int period() const;
    Real mass() const;
    Real exactMass() const;
    Real ionizationEnergy() const;
    Real electronAffinity() const;
    Real electronegativity() const;
    Real covalentRadius() const;
    Real vanDerWaalsRadius() const;
    Real boilingPoint() const;
    Real meltingPoint() const;
    int expectedValence() const;
    bool isValid() const;
    bool isMetal() const;
    bool isNonmetal() const;

    // operators
    inline bool operator==(const Element &element) const;
    inline bool operator!=(const Element &element) const;

    // static methods
    static Element fromName(const std::string &name);
    static Element fromName(const char *name);
    static Element fromSymbol(const std::string &symbol);
    static Element fromSymbol(const char *symbol);
    static Element fromSymbol(const char *symbol, size_t length);
    static Element fromSymbol(char symbol);
    static bool isValidAtomicNumber(AtomicNumberType atomicNumber);
    static bool isValidSymbol(const std::string &symbol);

private:
    AtomicNumberType m_atomicNumber;
};

} // end chemkit namespace

#include "element-inline.h"

#endif // CHEMKIT_ELEMENT_H
