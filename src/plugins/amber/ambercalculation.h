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

#ifndef AMBERCALCULATION_H
#define AMBERCALCULATION_H

#include <QtCore>

#include <chemkit/forcefieldcalculation.h>

class AmberParameters;

class AmberCalculation : public chemkit::ForceFieldCalculation
{
    public:
        virtual bool setup(const AmberParameters *parameters) = 0;

    protected:
        AmberCalculation(int type, int atomCount, int parameterCount);
};

class AmberBondCalculation : public AmberCalculation
{
    public:
        AmberBondCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup(const AmberParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class AmberAngleCalculation : public AmberCalculation
{
    public:
        AmberAngleCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c);

        bool setup(const AmberParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class AmberTorsionCalculation : public AmberCalculation
{
    public:
        AmberTorsionCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d);

        bool setup(const AmberParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

class AmberNonbondedCalculation : public AmberCalculation
{
    public:
        AmberNonbondedCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup(const AmberParameters *parameters);
        chemkit::Float energy() const;
        std::vector<chemkit::Vector3> gradient() const;
};

#endif // AMBERCALCULATION_H
