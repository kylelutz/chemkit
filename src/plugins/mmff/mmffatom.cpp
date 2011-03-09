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

#include "mmffatom.h"

#include "mmffforcefield.h"
#include "mmffparameters.h"

// --- Construction and Destruction ---------------------------------------- //
MmffAtom::MmffAtom(chemkit::ForceField *forceField, const chemkit::Atom *atom)
    : chemkit::ForceFieldAtom(forceField, atom)
{
    m_typeNumber = 0;
    m_formalCharge = 0;
}

// --- Properties ---------------------------------------------------------- //
const MmffForceField* MmffAtom::forceField() const
{
    return static_cast<const MmffForceField *>(chemkit::ForceFieldAtom::forceField());
}

void MmffAtom::setType(int typeNumber, chemkit::Float formalCharge)
{
    m_typeNumber = typeNumber;
    m_formalCharge = formalCharge;
}

QString MmffAtom::type() const
{
    return QString::number(m_typeNumber);
}

int MmffAtom::typeNumber() const
{
    return m_typeNumber;
}

chemkit::Float MmffAtom::formalCharge() const
{
    return m_formalCharge;
}

int MmffAtom::period() const
{
    return atom()->element().period();
}

// --- Parameters ---------------------------------------------------------- //
const MmffAtomParameters* MmffAtom::parameters() const
{
    const MmffParameters *mmffParameters = static_cast<const MmffForceField *>(forceField())->parameters();
    if(!mmffParameters)
        return 0;

    return mmffParameters->atomParameters(this);
}
