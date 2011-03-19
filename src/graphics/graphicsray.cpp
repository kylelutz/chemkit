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

#include "graphicsray.h"

namespace chemkit {

// === GraphicsRay ========================================================= //
/// \class GraphicsRay graphicsray.h chemkit/graphicsray.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsRay class represents a ray with an origin and
///        direction.

// --- Construction and Destruction ---------------------------------------- //
GraphicsRay::GraphicsRay()
    : m_origin(Point3f()),
      m_direction(-Vector3g::Z())
{
}

/// Creates a new ray with \p origin and \p direction.
GraphicsRay::GraphicsRay(const Point3f &origin, const Vector3g &direction)
    : m_origin(origin),
      m_direction(direction.normalized())
{
}

/// Creates a new ray with \p origin that points towards \p point.
GraphicsRay::GraphicsRay(const Point3f &origin, const Point3f &point)
    : m_origin(origin),
      m_direction((point - origin).normalized())
{
}

/// Destroys the graphics ray object.
GraphicsRay::~GraphicsRay()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the origin to \p origin.
void GraphicsRay::setOrigin(const Point3f &origin)
{
    m_origin = origin;
}

/// Returns the origin.
Point3f GraphicsRay::origin() const
{
    return m_origin;
}

/// Sets the direction to \p direction. The direction vector will be
/// normalized.
void GraphicsRay::setDirection(const Vector3g &direction)
{
    m_direction = direction.normalized();
}

/// Returns the direction vector.
Vector3g GraphicsRay::direction() const
{
    return m_direction;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the point at \p distance from the origin along the ray.
Point3f GraphicsRay::pointAt(float distance) const
{
    return m_origin.movedBy(distance, direction());
}

bool GraphicsRay::intersectsSphere(const Point3f &center, float radius, float *distance) const
{
    Vector3g dst = center - m_origin;

    float B = dst.dot(m_direction);
    float C = dst.dot(dst) - (radius * radius);
    float D = B*B - C;

    if(D < 0){
        return false;
    }

    if(distance){
        *distance = qMin(qAbs(-B + sqrt(D)), qAbs(-B - sqrt(D)));
    }

    return true;
}

bool GraphicsRay::intersectsCylinder(const Point3f &a, const Point3f &b, float radius, float *distance) const
{
    Vector3g ao = m_origin - a;
    Vector3g ab = b - a;
    Vector3g aoxab = ao.cross(ab);
    Vector3g vxab = direction().cross(ab);

    float A = vxab.dot(vxab);
    float B = 2 * vxab.dot(aoxab);
    float C = aoxab.dot(aoxab) - ab.dot(ab) * (radius * radius);
    float D = B*B - 4*A*C;

    if(D < 0){
        // no intersection
        return false;
    }

    float t = qMin((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));

    Point3f ip = pointAt(t);
    Point3f ip1 = ip - a;
    Point3f ip2 = ip - b;

    if(ip1.dot(ab) < 0 || ip2.dot(ab) > 0){
        // intersection below base or above top of the cylinder
        return false;
    }

    if(distance){
        *distance = t;
    }

    return true;
}

} // end chemkit namespace
