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

#include <chemkit/genericquaternion.h>

#include "point3g.h"
#include "vector3g.h"

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT Quaterniong : public GenericQuaternion<GraphicsFloat>
{
    public:
        // construction and destruction
        Quaterniong(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z, GraphicsFloat r);
        Quaterniong(const Point3g &point, GraphicsFloat r);
        Quaterniong(const Vector3g &vector, GraphicsFloat r);
        Quaterniong(const GenericQuaternion<GraphicsFloat> &quaternion);
        Quaterniong(const StaticVector<GraphicsFloat, 4> &quaternion);

        // properties
        Point3g toPoint3() const;
        Vector3g toVector3() const;

        // static methods
        static Quaterniong rotation(const Vector3g &axis, GraphicsFloat angle);
        static Quaterniong rotationRadians(const Vector3g &axis, GraphicsFloat angle);
        static Point3g rotate(const Point3g &point, const Vector3g &axis, GraphicsFloat angle);
        static Point3g rotateRadians(const Point3g &point, const Vector3g &axis, GraphicsFloat angle);
        static Vector3g rotate(const Vector3g &vector, const Vector3g &axis, GraphicsFloat angle);
        static Vector3g rotateRadians(const Vector3g &vector, const Vector3g &axis, GraphicsFloat angle);
};

} // end chemkit namespace

#include "quaterniong-inline.h"

#endif // CHEMKIT_QUATERNIONG_H
