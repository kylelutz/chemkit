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

#ifndef CHEMKIT_GRAPHICSRAY_H
#define CHEMKIT_GRAPHICSRAY_H

#include "graphics.h"

#include "vector3g.h"

#include <chemkit/point3.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsRay
{
    public:
        // construction and destruction
        GraphicsRay();
        GraphicsRay(const Point3f &origin, const Vector3g &direction);
        GraphicsRay(const Point3f &origin, const Point3f &point);
        ~GraphicsRay();

        // properties
        void setOrigin(const Point3f &origin);
        Point3f origin() const;
        void setDirection(const Vector3g &direction);
        Vector3g direction() const;

        // geometry
        Point3f pointAt(float distance) const;
        bool intersectsSphere(const Point3f &center, float radius, float *distance = 0) const;
        bool intersectsCylinder(const Point3f &a, const Point3f &b, float radius, float *distance = 0) const;

    private:
        Point3f m_origin;
        Vector3g m_direction;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSRAY_H
