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

#ifndef CHEMKIT_POINT3G_INLINE_H
#define CHEMKIT_POINT3G_INLINE_H

#include "point3g.h"

namespace chemkit {

// === Point3g ============================================================= //
/// \class Point3g point3g.h chemkit/point3g.h
/// \ingroup chemkit-graphics
/// \brief The Point3g class represent a location in
///        three-dimensional space.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new point at (\c 0, \c 0, \c 0).
inline Point3g::Point3g()
    : GenericPoint<GraphicsFloat>()
{
}

/// Create a new point at (\p x, \p y, \p z).
inline Point3g::Point3g(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
    : GenericPoint<GraphicsFloat>(x, y, z)
{
}

inline Point3g::Point3g(const GenericPoint<float> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

inline Point3g::Point3g(const GenericPoint<double> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

inline Point3g::Point3g(const StaticVector<float, 3> &point)
    : GenericPoint<GraphicsFloat>(point)
{

}

inline Point3g::Point3g(const StaticVector<double, 3> &point)
    : GenericPoint<GraphicsFloat>(point)
{
}

} // end chemkit namespace

#endif // CHEMKIT_POINT3G_INLINE_H
