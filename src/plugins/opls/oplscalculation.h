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

#ifndef OPLSCALCULATION_H
#define OPLSCALCULATION_H

#include <chemkit/forcefieldcalculation.h>

#include "oplsparameters.h"

class OplsCalculation : public chemkit::ForceFieldCalculation
{
public:
    virtual bool setup(const OplsParameters *parameters) = 0;

protected:
    OplsCalculation(int type, int atomCount, int parameterCount);
};

class OplsBondStrechCalculation : public OplsCalculation
{
public:
    OplsBondStrechCalculation(size_t a, size_t b);

    bool setup(const OplsParameters *parameters);
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class OplsAngleBendCalculation : public OplsCalculation
{
public:
    OplsAngleBendCalculation(size_t a, size_t b, size_t c);

    bool setup(const OplsParameters *parameters);
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class OplsTorsionCalculation : public OplsCalculation
{
public:
    OplsTorsionCalculation(size_t a, size_t b, size_t c, size_t d);

    bool setup(const OplsParameters *parameters);
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

class OplsNonbondedCalculation : public OplsCalculation
{
public:
    OplsNonbondedCalculation(size_t a, size_t b);

    bool setup(const OplsParameters *parameters);
    chemkit::Real energy(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
    std::vector<chemkit::Vector3> gradient(const chemkit::CartesianCoordinates *coordinates) const CHEMKIT_OVERRIDE;
};

#endif // OPLSCALCULATION_H
