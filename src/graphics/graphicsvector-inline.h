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

#ifndef CHEMKIT_GRAPHICSVECTOR_INLINE_H
#define CHEMKIT_GRAPHICSVECTOR_INLINE_H

#include "graphicsvector.h"

namespace chemkit {

// === GraphicsVector ====================================================== //
/// \class GraphicsVector graphicsvector.h chemkit/graphicsvector.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsVector class represents a direction in
///        three-dimensional space.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new vector containing (\c 0, \c 0, \c 0).
inline GraphicsVector::GraphicsVector()
    : GenericVector<GraphicsFloat>()
{
}

/// Create a new vector containing (\p x, \p y, \p z).
inline GraphicsVector::GraphicsVector(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
    : GenericVector<GraphicsFloat>(x, y, z)
{
}

inline GraphicsVector::GraphicsVector(const GenericVector<float> &vector)
    : GenericVector<GraphicsFloat>(vector)
{
}

inline GraphicsVector::GraphicsVector(const GenericVector<double> &vector)
    : GenericVector<GraphicsFloat>(vector)
{
}

inline GraphicsVector::GraphicsVector(const StaticVector<float, 3> &vector)
    : GenericVector<GraphicsFloat>(vector)
{
}

inline GraphicsVector::GraphicsVector(const StaticVector<double, 3> &vector)
    : GenericVector<GraphicsFloat>(vector)
{
}

// --- Properties ---------------------------------------------------------- //
inline Vector3 GraphicsVector::toVector3() const
{
    return Vector3(x(), y(), z());
}

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSVECTOR_INLINE_H
