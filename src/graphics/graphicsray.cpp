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
    : m_origin(GraphicsPoint()),
      m_direction(-GraphicsVector::Z())
{
}

/// Creates a new ray with \p origin and \p direction.
GraphicsRay::GraphicsRay(const GraphicsPoint &origin, const GraphicsVector &direction)
    : m_origin(origin),
      m_direction(direction.normalized())
{
}

/// Creates a new ray with \p origin that points towards \p point.
GraphicsRay::GraphicsRay(const GraphicsPoint &origin, const GraphicsPoint &point)
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
void GraphicsRay::setOrigin(const GraphicsPoint &origin)
{
    m_origin = origin;
}

/// Returns the origin.
GraphicsPoint GraphicsRay::origin() const
{
    return m_origin;
}

/// Sets the direction to \p direction. The direction vector will be
/// normalized.
void GraphicsRay::setDirection(const GraphicsVector &direction)
{
    m_direction = direction.normalized();
}

/// Returns the direction vector.
GraphicsVector GraphicsRay::direction() const
{
    return m_direction;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the point at \p distance from the origin along the ray.
GraphicsPoint GraphicsRay::pointAt(GraphicsFloat distance) const
{
    return m_origin.movedBy(distance, direction());
}

bool GraphicsRay::intersectsSphere(const GraphicsPoint &center, GraphicsFloat radius, GraphicsFloat *distance) const
{
    GraphicsVector dst = center - m_origin;

    GraphicsFloat B = dst.dot(m_direction);
    GraphicsFloat C = dst.dot(dst) - (radius * radius);
    GraphicsFloat D = B*B - C;

    if(D < 0){
        return false;
    }

    if(distance){
        *distance = qMin(qAbs(-B + sqrt(D)), qAbs(-B - sqrt(D)));
    }

    return true;
}

bool GraphicsRay::intersectsCylinder(const GraphicsPoint &a, const GraphicsPoint &b, GraphicsFloat radius, GraphicsFloat *distance) const
{
    GraphicsVector ao = m_origin - a;
    GraphicsVector ab = b - a;
    GraphicsVector aoxab = ao.cross(ab);
    GraphicsVector vxab = direction().cross(ab);

    GraphicsFloat A = vxab.dot(vxab);
    GraphicsFloat B = 2 * vxab.dot(aoxab);
    GraphicsFloat C = aoxab.dot(aoxab) - ab.dot(ab) * (radius * radius);
    GraphicsFloat D = B*B - 4*A*C;

    if(D < 0){
        // no intersection
        return false;
    }

    GraphicsFloat t = qMin((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));

    GraphicsPoint ip = pointAt(t);
    GraphicsPoint ip1 = ip - a;
    GraphicsPoint ip2 = ip - b;

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
