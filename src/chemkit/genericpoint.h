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

#ifndef CHEMKIT_GENERICPOINT_H
#define CHEMKIT_GENERICPOINT_H

#include "staticvector.h"

namespace chemkit {

template<typename T>
class GenericPoint : public StaticVector<T, 3>
{
    public:
        // construction and destruction
        GenericPoint();
        GenericPoint(T x, T y, T z);
        GenericPoint(const StaticVector<float, 3> &vector);
        GenericPoint(const StaticVector<double, 3> &vector);

        // properties
        T x() const;
        T y() const;
        T z() const;
        void moveBy(T dx, T dy, T dz);
        void moveBy(const StaticVector<T, 3> &vector);
        void moveBy(T distance, const StaticVector<T, 3> &vector);
        GenericPoint<T> movedBy(T dx, T dy, T dz) const;
        GenericPoint<T> movedBy(const StaticVector<T, 3> &vector) const;
        GenericPoint<T> movedBy(T distance, const StaticVector<T, 3> &direction) const;

        // math
        T distance(const GenericPoint<T> &point) const;
        GenericPoint<T> midpoint(const GenericPoint<T> &point) const;

        // static methods
        static T distance(const GenericPoint<T> &a, const GenericPoint<T> &b);
        static T distanceSquared(const GenericPoint<T> &a, const GenericPoint<T> &b);
        static T angle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c);
        static T angleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c);
        static T torsionAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T torsionAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T wilsonAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T wilsonAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static GenericPoint<T> midpoint(const GenericPoint<T> &a, const GenericPoint<T> &b);
};

} // end chemkit namespace

#include "genericpoint-inline.h"

#endif // CHEMKIT_GENERICPOINT_H
