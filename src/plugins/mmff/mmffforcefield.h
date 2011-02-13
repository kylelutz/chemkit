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

#ifndef MMFFFORCEFIELD_H
#define MMFFFORCEFIELD_H

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>

#include "mmffcalculation.h"

class MmffAtom;
class MmffParameters;

class MmffForceField : public chemkit::ForceField
{
    public:
        // construction and destruction
        MmffForceField();
        ~MmffForceField();

        // atoms
        MmffAtom* atom(const chemkit::Atom *atom);
        const MmffAtom* atom(const chemkit::Atom *atom) const;

        // parameterization
        virtual bool setup();
        const MmffParameters* parameters() const;

        // static methods
        static bool isAromatic(const chemkit::Ring *ring);
        static bool isAromatic(const chemkit::Atom *atom);
        static bool isAromatic(const chemkit::Bond *bond);
        static int piElectronCount(const chemkit::Ring *ring);

    private:
        int ringPosition(const chemkit::Atom *atom, const chemkit::Ring *ring) const;
        bool atomsWithinTwoBonds(const chemkit::Atom *a, const chemkit::Atom *b);
        bool assignPartialCharge(MmffAtom *a);
        static void clearParametersCache();

    private:
        MmffParameters *m_parameters;

        static QHash<QString, MmffParameters *> m_parametersCache;
};

#endif // MMFFFORCEFIELD_H
