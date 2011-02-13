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

#ifndef CHEMKIT_GRAPHICSPOINT_INLINE_H
#define CHEMKIT_GRAPHICSPOINT_INLINE_H

#include "graphicspoint.h"

namespace chemkit {

// === GraphicsPoint ======================================================= //
/// \class GraphicsPoint graphicspoint.h chemkit/graphicspoint.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsPoint class represent a location in
///        three-dimensional space.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new point at (\c 0, \c 0, \c 0).
inline GraphicsPoint::GraphicsPoint()
    : GenericPoint<GraphicsFloat>()
{
}

/// Create a new point at (\p x, \p y, \p z).
inline GraphicsPoint::GraphicsPoint(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
    : GenericPoint<GraphicsFloat>(x, y, z)
{
}

inline GraphicsPoint::GraphicsPoint(const GenericPoint<float> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

inline GraphicsPoint::GraphicsPoint(const GenericPoint<double> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

inline GraphicsPoint::GraphicsPoint(const StaticVector<float, 3> &point)
    : GenericPoint<GraphicsFloat>(point)
{

}

inline GraphicsPoint::GraphicsPoint(const StaticVector<double, 3> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

// --- Properties ---------------------------------------------------------- //
inline Point GraphicsPoint::toPoint() const
{
    return Point(x(), y(), z());
}

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPOINT_INLINE_H
