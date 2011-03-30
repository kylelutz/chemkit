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

#ifndef MMFFATOM_H
#define MMFFATOM_H

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>
#include <chemkit/forcefieldatom.h>

class MmffForceField;
struct MmffAtomParameters;

class MmffAtom : public chemkit::ForceFieldAtom
{
    public:
        // construction and destruction
        MmffAtom(chemkit::ForceField *forceField, const chemkit::Atom *atom);

        // properties
        const MmffForceField* forceField() const;
        void setType(int typeNumber, chemkit::Float formalCharge = 0);
        virtual std::string type() const;
        int typeNumber() const;
        chemkit::Float formalCharge() const;
        int period() const;

        // parameters
        const MmffAtomParameters* parameters() const;

    private:
        int m_typeNumber;
        chemkit::Float m_formalCharge;
};

#endif // MMFFATOM_H
