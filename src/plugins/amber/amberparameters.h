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

#ifndef AMBERPARAMETERS_H
#define AMBERPARAMETERS_H

#include <chemkit/forcefieldatom.h>

struct AmberBondParameters
{
    chemkit::Float kb;
    chemkit::Float r0;
};

struct AmberAngleParameters
{
    chemkit::Float ka;
    chemkit::Float theta0;
};

struct AmberTorsionParameters
{
    chemkit::Float V1;
    chemkit::Float V2;
    chemkit::Float V3;
    chemkit::Float V4;
    chemkit::Float gamma1;
    chemkit::Float gamma2;
    chemkit::Float gamma3;
    chemkit::Float gamma4;
};

struct AmberNonbondedParameters
{
    chemkit::Float vanDerWaalsRadius;
    chemkit::Float wellDepth;
};

class AmberParameters
{
    public:
        // construction and destruction
        AmberParameters();

        // parameters
        const AmberBondParameters* bondParameters(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b) const;
        const AmberAngleParameters* angleParameters(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c) const;
        const AmberTorsionParameters* torsionParameters(const chemkit::ForceFieldAtom *a, const chemkit::ForceFieldAtom *b, const chemkit::ForceFieldAtom *c, const chemkit::ForceFieldAtom *d) const;
        const AmberNonbondedParameters* nonbondedParameters(const chemkit::ForceFieldAtom *atom) const;
};

#endif // AMBERPARAMETERS_H
