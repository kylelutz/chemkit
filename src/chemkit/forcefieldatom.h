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

#ifndef CHEMKIT_FORCEFIELDATOM_H
#define CHEMKIT_FORCEFIELDATOM_H

#include "chemkit.h"

#include "point3.h"
#include "vector3.h"

namespace chemkit {

class Atom;
class ForceField;
class ForceFieldAtomPrivate;

class CHEMKIT_EXPORT ForceFieldAtom
{
public:
    // construction and destruction
    ForceFieldAtom(ForceField *forceField, const Atom *atom);
    virtual ~ForceFieldAtom();

    // properties
    const Atom* atom() const;
    int index() const;
    virtual bool setType(const std::string &type);
    virtual std::string type() const;
    void setCharge(Real charge);
    Real charge() const;
    bool isSetup() const;
    ForceField* forceField() const;

    // calculations
    Real energy() const;
    Vector3 gradient() const;

    // structure
    bool isOneFour(const ForceFieldAtom *atom) const;

    // geometry
    void setPosition(const Point3 &position);
    Point3 position() const;
    void moveBy(const Vector3 &vector);
    void moveBy(Real dx, Real dy, Real dz);

private:
    ForceFieldAtomPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDATOM_H
