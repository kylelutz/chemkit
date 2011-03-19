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

#ifndef CHEMKIT_QUATERNIONG_H
#define CHEMKIT_QUATERNIONG_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/genericquaternion.h>

#include "vector3g.h"

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT Quaterniong : public GenericQuaternion<float>
{
    public:
        // construction and destruction
        Quaterniong(float x, float y, float z, float r);
        Quaterniong(const Point3f &point, float r);
        Quaterniong(const Vector3g &vector, float r);
        Quaterniong(const GenericQuaternion<float> &quaternion);
        Quaterniong(const StaticVector<float, 4> &quaternion);

        // properties
        Point3f toPoint3() const;
        Vector3g toVector3() const;

        // static methods
        static Quaterniong rotation(const Vector3g &axis, float angle);
        static Quaterniong rotationRadians(const Vector3g &axis, float angle);
        static Point3f rotate(const Point3f &point, const Vector3g &axis, float angle);
        static Point3f rotateRadians(const Point3f &point, const Vector3g &axis, float angle);
        static Vector3g rotate(const Vector3g &vector, const Vector3g &axis, float angle);
        static Vector3g rotateRadians(const Vector3g &vector, const Vector3g &axis, float angle);
};

} // end chemkit namespace

#include "quaterniong-inline.h"

#endif // CHEMKIT_QUATERNIONG_H
