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

#ifndef CHEMKIT_FORCEFIELDATOM_H
#define CHEMKIT_FORCEFIELDATOM_H

#include "chemkit.h"

#include <QtCore>

#include "point.h"
#include "vector.h"

namespace chemkit {

class Atom;
class ForceField;
class ForceFieldAtomPrivate;

class CHEMKIT_EXPORT ForceFieldAtom
{
    public:
        // construction and destruction
        ForceFieldAtom(ForceField *forceField, const Atom *atom);
        virtual ~ForceFieldAtom();

        // properties
        const Atom* atom() const;
        int index() const;
        virtual bool setType(const QString &type);
        virtual QString type() const;
        void setCharge(Float charge);
        Float charge() const;
        bool isSetup() const;
        ForceField* forceField();
        const ForceField* forceField() const;

        // calculations
        Float energy() const;
        Vector gradient() const;

        // structure
        bool isOneFour(const ForceFieldAtom *atom) const;

        // geometry
        void setPosition(const Point &position);
        Point position() const;
        void moveBy(const Vector &vector);
        void moveBy(Float dx, Float dy, Float dz);

    private:
        ForceFieldAtomPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDATOM_H
