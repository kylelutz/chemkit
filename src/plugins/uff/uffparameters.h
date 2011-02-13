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

#ifndef UFFPARAMETERS_H
#define UFFPARAMETERS_H

#include <chemkit/forcefieldatom.h>

struct UffAtomParameters {
    const char *type;
    chemkit::Float r; // bond length (angstroms)
    chemkit::Float theta; // angle (degrees)
    chemkit::Float x; // distance (angstroms)
    chemkit::Float D; // energy (kcal/mol)
    chemkit::Float zeta; // scale
    chemkit::Float Z; // charge
    chemkit::Float V; // torsional barrier (kcal/mol)
    chemkit::Float U;
    chemkit::Float X; // electronegativity (eV)
    chemkit::Float hard;
    chemkit::Float radius;
};

class UffParameters
{
    public:
        // construction and destruction
        UffParameters();
        ~UffParameters();

        // parameters
        const UffAtomParameters* parameters(const chemkit::ForceFieldAtom *atom) const;
};

#endif // UFFPARAMETERS_H
