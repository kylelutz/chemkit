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

#ifndef CHEMKIT_ATOMTYPER_H
#define CHEMKIT_ATOMTYPER_H

#include "chemkit.h"

#include <string>
#include <vector>

#include "variant.h"

namespace chemkit {

class Atom;
class Molecule;
class AtomTyperPrivate;

class CHEMKIT_EXPORT AtomTyper
{
public:
    // typedefs
    typedef AtomTyper* (*CreateFunction)();

    // construction and destruction
    virtual ~AtomTyper();

    // properties
    std::string name() const;
    void setMolecule(const Molecule *molecule);
    const Molecule* molecule() const;

    // types
    virtual Variant type(int index) const;
    virtual Variant type(const Atom *atom) const;
    virtual int typeNumber(int index) const;
    virtual int typeNumber(const Atom *atom) const;
    virtual std::string typeString(int index) const;
    virtual std::string typeString(const Atom *atom) const;

    // static methods
    static AtomTyper* create(const std::string &name);
    static std::vector<std::string> typers();

protected:
    AtomTyper(const std::string &name);
    virtual void assignTypes(const Molecule *molecule);

private:
    AtomTyperPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_ATOMTYPER_H
