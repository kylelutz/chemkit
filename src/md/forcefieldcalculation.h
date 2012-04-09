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

#include <boost/shared_ptr.hpp>

#include <chemkit/vector3.h>

namespace chemkit {

class Topology;
class ForceField;
class CartesianCoordinates;
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
    boost::shared_ptr<Topology> topology() const;

    // atoms
    size_t atom(size_t index) const;
    std::vector<size_t> atoms() const;
    size_t atomCount() const;
    std::string atomType(size_t index) const;

    // parameters
    void setParameter(int index, Real value);
    Real parameter(int index) const;
    std::vector<Real> parameters() const;
    int parameterCount() const;

    // calculations
    virtual Real energy(const CartesianCoordinates *coordinates) const;
    virtual std::vector<Vector3> gradient(const CartesianCoordinates *coordinates) const;
    std::vector<Vector3> numericalGradient(const CartesianCoordinates *coordinates) const;

protected:
    ForceFieldCalculation(int type, size_t atomCount, size_t parameterCount);
    virtual ~ForceFieldCalculation();
    void setAtom(size_t index, size_t atom);

private:
    void setSetup(bool setup);
    void setForceField(ForceField *forceField);

    friend class ForceField;

private:
    ForceFieldCalculationPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDCALCULATION_H
