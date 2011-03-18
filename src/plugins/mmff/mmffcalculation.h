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

#ifndef MMFFCALCULATION_H
#define MMFFCALCULATION_H

#include <chemkit/forcefieldcalculation.h>

class MmffAtom;
class MmffParameters;

class MmffCalculation : public chemkit::ForceFieldCalculation
{
    public:
        const MmffAtom* atom(int index) const;
        virtual bool setup(const MmffParameters *parameters) = 0;

    protected:
        MmffCalculation(int type, int atomCount, int parameterCount);
};

class MmffBondStrechCalculation : public MmffCalculation
{
    public:
        MmffBondStrechCalculation(const MmffAtom *a, const MmffAtom *b);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffAngleBendCalculation : public MmffCalculation
{
    public:
        MmffAngleBendCalculation(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffStrechBendCalculation : public MmffCalculation
{
    public:
        MmffStrechBendCalculation(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffOutOfPlaneBendingCalculation : public MmffCalculation
{
    public:
        MmffOutOfPlaneBendingCalculation(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffTorsionCalculation : public MmffCalculation
{
    public:
        MmffTorsionCalculation(const MmffAtom *a, const MmffAtom *b, const MmffAtom *c, const MmffAtom *d);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffVanDerWaalsCalculation : public MmffCalculation
{
    public:
        MmffVanDerWaalsCalculation(const MmffAtom *a, const MmffAtom *b);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

class MmffElectrostaticCalculation : public MmffCalculation
{
    public:
        MmffElectrostaticCalculation(const MmffAtom *a, const MmffAtom *b);

        bool setup(const MmffParameters *parameters);
        chemkit::Float energy() const;
        QVector<chemkit::Vector3> gradient() const;
};

#endif // MMFFCALCULATION_H
