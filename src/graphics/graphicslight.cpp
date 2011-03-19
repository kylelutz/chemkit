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

#include "graphicslight.h"

namespace chemkit {

// === GraphicsLightPrivate ================================================ //
class GraphicsLightPrivate
{
    public:
        Point3f position;
        Vector3f direction;
};

// === GraphicsLight ======================================================= //
/// \class GraphicsLight graphicslight.h chemkit/graphicslight.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsLight class represents a light in a graphics
///        view.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics light.
GraphicsLight::GraphicsLight()
    : d(new GraphicsLightPrivate)
{
}

/// Destroys the graphics light object.
GraphicsLight::~GraphicsLight()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the position of the light to \p position.
void GraphicsLight::setPosition(const Point3f &position)
{
    d->position = position;
}

/// Returns the position of the light.
Point3f GraphicsLight::position() const
{
    return d->position;
}

/// Sets the direction the camera points to \p direction.
void GraphicsLight::setDirection(const Vector3f &direction)
{
    d->direction = direction;
}

/// Returns the direction that the camera points.
Vector3f GraphicsLight::direction() const
{
    return d->direction;
}

} // end chemkit namespace
