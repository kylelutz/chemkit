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

#ifndef UFFCALCULATION_H
#define UFFCALCULATION_H

#include <chemkit/forcefieldcalculation.h>

#include "uffparameters.h"

class UffCalculation : public chemkit::ForceFieldCalculation
{
public:
    UffCalculation(int type, int atomCount, int parameterCount);

    virtual bool setup() = 0;

protected:
    chemkit::Real bondOrder(size_t a,  size_t b) const;
    chemkit::Real bondLength(const UffAtomParameters *a, const UffAtomParameters *b, chemkit::Real bondOrder) const;
    const UffAtomParameters* parameters(const std::string &type) const;
};

class UffBondStrechCalculation : public UffCalculation
{
public:
    UffBondStrechCalculation(size_t a, size_t b);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class UffAngleBendCalculation : public UffCalculation
{
public:
    UffAngleBendCalculation(size_t a, size_t b, size_t c);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class UffTorsionCalculation : public UffCalculation
{
public:
    UffTorsionCalculation(size_t a, size_t b, size_t c, size_t d);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class UffInversionCalculation : public UffCalculation
{
public:
    UffInversionCalculation(size_t a, size_t b, size_t c, size_t d);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class UffVanDerWaalsCalculation : public UffCalculation
{
public:
    UffVanDerWaalsCalculation(size_t a, size_t b);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class UffElectrostaticCalculation : public UffCalculation
{
public:
    UffElectrostaticCalculation(size_t a, size_t b);

    bool setup();
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

#endif // UFFCALCULATION_H
