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

#ifndef MMFFATOMTYPER_H
#define MMFFATOMTYPER_H

#include <chemkit/ring.h>
#include <chemkit/atomtyper.h>

class MmffAtomTyper : public chemkit::AtomTyper
{
    public:
        // construction and destruction
        MmffAtomTyper(const chemkit::Molecule *molecule = 0);
        ~MmffAtomTyper();

        // types
        int typeNumber(int index) const;
        int typeNumber(const chemkit::Atom *atom) const;

        // charges
        chemkit::Float formalCharge(int index) const;
        chemkit::Float formalCharge(const chemkit::Atom *atom) const;

    protected:
        void assignTypes(const chemkit::Molecule *molecule);

    private:
        void setType(int index, int type, chemkit::Float formalCharge = 0);
        void setType(int index, const chemkit::Atom *atom);
        void setHydrogenType(int index, const chemkit::Atom *atom);
        void setCarbonType(int index, const chemkit::Atom *atom);
        void setNitrogenType(int index, const chemkit::Atom *atom);
        void setOxygenType(int index, const chemkit::Atom *atom);
        void setSulfurType(int index, const chemkit::Atom *atom);
        void setAromaticType(int index, const chemkit::Atom *atom, const chemkit::Ring *ring, int position);

    private:
        QVector<int> m_types;
        QVector<chemkit::Float> m_formalCharges;
        QVector<chemkit::Float> m_partialCharges;
};

#endif // MMFFATOMTYPER_H
