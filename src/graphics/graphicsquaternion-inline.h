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

#ifndef CHEMKIT_GRAPHICSQUATERNION_INLINE_H
#define CHEMKIT_GRAPHICSQUATERNION_INLINE_H

#include "graphicsquaternion.h"

#include <chemkit/constants.h>

namespace chemkit {

// === GraphicsQuaterion =================================================== //
/// \class GraphicsQuaternion graphicsquaternion.h chemkit/graphicsquaternion.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsQuaternion class represents a quaternion.

// --- Construction and Destruction ---------------------------------------- //    
inline GraphicsQuaternion::GraphicsQuaternion(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(x, y, z, r)
{
}

inline GraphicsQuaternion::GraphicsQuaternion(const Point3g &point, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(point.x(), point.y(), point.z(), r)
{
}

inline GraphicsQuaternion::GraphicsQuaternion(const GraphicsVector &vector, GraphicsFloat r)
    : GenericQuaternion<GraphicsFloat>(vector.x(), vector.y(), vector.z(), r)
{
}

inline GraphicsQuaternion::GraphicsQuaternion(const GenericQuaternion<GraphicsFloat> &quaternion)
    : GenericQuaternion<GraphicsFloat>(quaternion)
{
}

inline GraphicsQuaternion::GraphicsQuaternion(const StaticVector<GraphicsFloat, 4> &quaternion)
    : GenericQuaternion<GraphicsFloat>(quaternion)
{
}

// --- Properties ---------------------------------------------------------- //
inline Point3g GraphicsQuaternion::toPoint3() const
{
    return Point3g(x(), y(), z());
}

inline GraphicsVector GraphicsQuaternion::toVector() const
{
    return GraphicsVector(x(), y(), z());
}

// --- Static Methods ------------------------------------------------------ //
inline GraphicsQuaternion GraphicsQuaternion::rotation(const GraphicsVector &axis, GraphicsFloat angle)
{
    return rotationRadians(axis, angle * chemkit::constants::DegreesToRadians);
}

inline GraphicsQuaternion GraphicsQuaternion::rotationRadians(const GraphicsVector &axis, GraphicsFloat angle)
{
    return GraphicsQuaternion(axis.x() * sin(angle/2.0),
                              axis.y() * sin(angle/2.0),
                              axis.z() * sin(angle/2.0),
                              cos(angle/2.0));
}

inline Point3g GraphicsQuaternion::rotate(const Point3g &point, const GraphicsVector &axis, GraphicsFloat angle)
{
    return rotateRadians(point, axis, angle * constants::DegreesToRadians);
}

inline Point3g GraphicsQuaternion::rotateRadians(const Point3g &point, const GraphicsVector &axis, GraphicsFloat angle)
{
    GraphicsQuaternion p(point.x(), point.y(), point.z(), 0);
    GraphicsQuaternion q = rotationRadians(axis, angle);
    GraphicsQuaternion qc = q.conjugate();

    GraphicsQuaternion qp = q.multiply(p);
    GraphicsQuaternion r = qp.multiply(qc);

    return r.toPoint3();
}

inline GraphicsVector GraphicsQuaternion::rotate(const GraphicsVector &vector, const GraphicsVector &axis, GraphicsFloat angle)
{
    return rotateRadians(vector, axis, angle * constants::DegreesToRadians);
}

inline GraphicsVector GraphicsQuaternion::rotateRadians(const GraphicsVector &vector, const GraphicsVector &axis, GraphicsFloat angle)
{
    GraphicsQuaternion p(vector.x(), vector.y(), vector.z(), 0);
    GraphicsQuaternion q = rotationRadians(axis, angle);
    GraphicsQuaternion qc = q.conjugate();

    GraphicsQuaternion qp = q.multiply(p);
    GraphicsQuaternion r = qp.multiply(qc);

    return r.toVector();
}

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSQUATERNION_INLINE_H
