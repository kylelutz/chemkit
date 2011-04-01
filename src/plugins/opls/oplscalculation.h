/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
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
        OplsBondStrechCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup(const OplsParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class OplsAngleBendCalculation : public OplsCalculation
{
    public:
        OplsAngleBendCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c);

        bool setup(const OplsParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class OplsTorsionCalculation : public OplsCalculation
{
    public:
        OplsTorsionCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d);

        bool setup(const OplsParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class OplsNonbondedCalculation : public OplsCalculation
{
    public:
        OplsNonbondedCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup(const OplsParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

#endif // OPLSCALCULATION_H
