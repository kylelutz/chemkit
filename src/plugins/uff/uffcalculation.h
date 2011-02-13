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
        chemkit::Float bondOrder(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b) const;
        chemkit::Float bondLength(const UffAtomParameters *a, const UffAtomParameters *b, chemkit::Float bondOrder) const;
        const UffAtomParameters* parameters(const chemkit::ForceFieldAtom *atom) const;
};

class UffBondStrechCalculation : public UffCalculation
{
    public:
        UffBondStrechCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup();
        chemkit::Float energy() const;
        QVector<chemkit::Vector> gradient() const;
};

class UffAngleBendCalculation : public UffCalculation
{
    public:
        UffAngleBendCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c);

        bool setup();
        chemkit::Float energy() const;
        QVector<chemkit::Vector> gradient() const;
};

class UffTorsionCalculation : public UffCalculation
{
    public:
        UffTorsionCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d);

        bool setup();
        chemkit::Float energy() const;
        QVector<chemkit::Vector> gradient() const;
};

class UffInversionCalculation : public UffCalculation
{
    public:
        UffInversionCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d);

        bool setup();
        chemkit::Float energy() const;
        QVector<chemkit::Vector> gradient() const;
};

class UffVanDerWaalsCalculation : public UffCalculation
{
    public:
        UffVanDerWaalsCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup();
        chemkit::Float energy() const;
        QVector<chemkit::Vector> gradient() const;
};

class UffElectrostaticCalculation : public UffCalculation
{
    public:
        UffElectrostaticCalculation(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b);

        bool setup();
        chemkit::Float energy() const;
};

#endif // UFFCALCULATION_H
