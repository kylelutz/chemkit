/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      m_direction(-Vector3f::UnitZ())
{
}

/// Creates a new ray with \p origin and \p direction.
GraphicsRay::GraphicsRay(const Point3f &origin, const Vector3f &direction)
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
void GraphicsRay::setDirection(const Vector3f &direction)
{
    m_direction = direction.normalized();
}

/// Returns the direction vector.
Vector3f GraphicsRay::direction() const
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
    Vector3f dst = center - m_origin;

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
    Vector3f ao = m_origin - a;
    Vector3f ab = b - a;
    Vector3f aoxab = ao.cross(ab);
    Vector3f vxab = direction().cross(ab);

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
