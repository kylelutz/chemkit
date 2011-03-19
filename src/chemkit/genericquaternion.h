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

#ifndef CHEMKIT_GENERICQUATERNION_H
#define CHEMKIT_GENERICQUATERNION_H

#include "chemkit.h"

#include "genericpoint.h"
#include "genericvector.h"
#include "staticvector.h"

namespace chemkit {

template<typename T>
class GenericQuaternion : public StaticVector<T, 4>
{
    public:
        // construction and destruction
        GenericQuaternion();
        GenericQuaternion(T x, T y, T z, T r);
        GenericQuaternion(const StaticVector<T, 4> &quaternion);

        // properties
        T x() const;
        T y() const;
        T z() const;
        T r() const;
        GenericPoint<T> toPoint3() const;
        GenericVector<T> toVector3() const;

        // math
        GenericQuaternion<T> multiply(const GenericQuaternion<T> &quaternion) const;
        GenericQuaternion<T> conjugate() const;

        // operators
        GenericQuaternion<T> operator*(const GenericQuaternion<T> &quaternion) const;

    private:
        static GenericQuaternion<T> multiply(const GenericQuaternion<T> &a, const GenericQuaternion<T> &b);
};

} // end chemkit namespace

#include "genericquaternion-inline.h"

#endif // CHEMKIT_GENERICQUATERNION_H
