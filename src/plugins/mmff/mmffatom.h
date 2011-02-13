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
        void setType();
        void setType(int typeNumber, chemkit::Float formalCharge = 0);
        void setHydrogenType(const MmffAtom *neighborAtom);
        void setCarbonType();
        void setNitrogenType();
        void setOxygenType();
        void setSulfurType();
        void setAromaticType(const chemkit::Ring *aromaticRing, int position);
        virtual QString type() const;
        int typeNumber() const;
        chemkit::Float formalCharge() const;
        int period() const;
        bool setCharge();

        // parameters
        const MmffAtomParameters* parameters() const;

    private:
        bool isGuanidinium(const chemkit::Atom *atom) const;
        bool isResonant(const chemkit::Atom *atom) const;
        bool isAmide(const chemkit::Atom *atom) const;
        bool isPhosphate(const chemkit::Atom *atom) const;
        bool isSulfate(const chemkit::Atom *atom) const;
        bool isThiocarboxylate(const chemkit::Atom *atom) const;
        bool isPositiveAromaticNitrogenRing(const chemkit::Ring *ring) const;

    private:
        int m_typeNumber;
        chemkit::Float m_formalCharge;
};

#endif // MMFFATOM_H
