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

#ifndef CHEMKIT_FORCEFIELDCALCULATION_H
#define CHEMKIT_FORCEFIELDCALCULATION_H

#include "md.h"

#include <vector>

#include <boost/array.hpp>

#include <chemkit/point3.h>
#include <chemkit/vector3.h>

namespace chemkit {

class ForceField;
class ForceFieldAtom;
class ForceFieldCalculationPrivate;

class CHEMKIT_MD_EXPORT ForceFieldCalculation
{
public:
    // enumerations
    enum Type {
        BondStrech = 0x01,
        AngleBend = 0x02,
        Torsion = 0x04,
        Inversion = 0x08,
        VanDerWaals = 0x10,
        Electrostatic = 0x20
    };

    // properties
    int type() const;
    bool isSetup() const;
    ForceField* forceField() const;

    // atoms
    const ForceFieldAtom* atom(int index) const;
    std::vector<const ForceFieldAtom *> atoms() const;
    int atomCount() const;
    bool contains(const ForceFieldAtom *atom) const;

    // parameters
    void setParameter(int index, Real value);
    Real parameter(int index) const;
    std::vector<Real> parameters() const;
    int parameterCount() const;

    // calculations
    virtual Real energy() const;
    virtual std::vector<Vector3> gradient() const;
    std::vector<Vector3> numericalGradient() const;

protected:
    ForceFieldCalculation(int type, int atomCount, int parameterCount);
    virtual ~ForceFieldCalculation();
    void setAtom(int index, const ForceFieldAtom *atom);
    inline Real distance(const ForceFieldAtom *a, const ForceFieldAtom *b) const;
    inline boost::array<Vector3, 2> distanceGradient(const ForceFieldAtom *a, const ForceFieldAtom *b) const;
    inline Real bondAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const;
    inline Real bondAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const;
    inline boost::array<Vector3, 3> bondAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const;
    inline boost::array<Vector3, 3> bondAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const;
    inline Real torsionAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline Real torsionAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline boost::array<Vector3, 4> torsionAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline boost::array<Vector3, 4> torsionAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline Real wilsonAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline Real wilsonAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline boost::array<Vector3, 4> wilsonAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;
    inline boost::array<Vector3, 4> wilsonAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const;

private:
    void setSetup(bool setup);

    friend class ForceField;

private:
    ForceFieldCalculationPrivate* const d;
};

} // end chemkit namespace

#include "forcefieldcalculation-inline.h"

#endif // CHEMKIT_FORCEFIELDCALCULATION_H
