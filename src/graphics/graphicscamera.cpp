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

#include "graphicscamera.h"

#include <chemkit/geometry.h>

namespace chemkit {

// === GraphicsCameraPrivate =============================================== //
class GraphicsCameraPrivate
{
public:
    Point3f position;
    Vector3f direction;
    Vector3f upVector;
    Point3f focus;
    bool changed;
};

// === GraphicsCamera ====================================================== //
/// \class GraphicsCamera graphicscamera.h chemkit/graphicscamera.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsCamera class represents a camera in a graphics
///        view.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics camera object.
GraphicsCamera::GraphicsCamera()
    : d(new GraphicsCameraPrivate)
{
    d->position = Point3f(0, 0, 0);
    d->direction = -Vector3f::UnitZ();
    d->upVector = Vector3f::UnitY();
    d->focus = Point3f(0, 0, 0);
    d->changed = true;
}

/// Creates a new graphics camera object at \p position.
GraphicsCamera::GraphicsCamera(const Point3f &position)
    : d(new GraphicsCameraPrivate)
{
    d->position = position;
    d->direction = -Vector3f::UnitZ();
    d->upVector = Vector3f::UnitY();
    d->focus = Point3f(0, 0, 0);
    d->changed = true;
}

/// \overload
GraphicsCamera::GraphicsCamera(float x, float y, float z)
    : d(new GraphicsCameraPrivate)
{
    d->position = Point3f(x, y, z);
    d->direction = -Vector3f::UnitZ();
    d->upVector = Vector3f::UnitY();
    d->focus = Point3f(0, 0, 0);
    d->changed = true;
}

/// Destroys the camera object.
GraphicsCamera::~GraphicsCamera()
{
    delete d;
}

// --- Position ------------------------------------------------------------ //
/// Sets the position of the camera to \p position.
void GraphicsCamera::setPosition(const Point3f &position)
{
    d->position = position;
    setChanged(true);
}

/// Sets the position of the camera to (\p x, \p y, \p z).
void GraphicsCamera::setPosition(float x, float y, float z)
{
    setPosition(Point3f(x, y, z));
}

/// Returns the position of the camera.
Point3f GraphicsCamera::position() const
{
    return d->position;
}

/// Returns the x component of the camera's position.
float GraphicsCamera::x() const
{
    return position().x();
}

/// Returns the y component of the camera's position.
float GraphicsCamera::y() const
{
    return position().y();
}

/// Returns the z component of the camera's position.
float GraphicsCamera::z() const
{
    return position().z();
}

/// Sets the camera's position to \p position.
void GraphicsCamera::moveTo(const Point3f &position)
{
    setPosition(position);
}

/// Sets the camera's position to (\p x, \p y, \p z).
void GraphicsCamera::moveTo(float x, float y, float z)
{
    setPosition(x, y, z);
}

/// Moves the camera by \p vector.
void GraphicsCamera::moveBy(const Vector3f &vector)
{
    setPosition(position() + vector);
}

/// Moves the camera by (\p dx, \p dy, \p dz).
void GraphicsCamera::moveBy(float dx, float dy, float dz)
{
    moveBy(Vector3f(dx, dy, dz));
}

/// Moves the camera by \p distance in \p direction.
void GraphicsCamera::moveBy(float distance, const Vector3f &direction)
{
    setPosition(position() + (direction.normalized() * distance));
}

/// Moves the camera fowards by \p distance.
void GraphicsCamera::moveFoward(float distance)
{
    moveBy(distance, direction());
}

/// Moves the camera backwards by \p distance.
void GraphicsCamera::moveBackward(float distance)
{
    moveBy(-distance, direction());
}

/// Rotates the camera around \p axis by \p angle degrees. If
/// \p rotateDirection is \c true the camera's direction will also
/// be rotated.
void GraphicsCamera::rotate(const Vector3f &axis, float angle, bool rotateDirection)
{
    setPosition(chemkit::geometry::rotate(position(), axis, angle));

    if(rotateDirection){
        setDirection(chemkit::geometry::rotate(direction(), axis, angle));
        setUpVector(chemkit::geometry::rotate(upVector(), axis, angle));
    }
}

/// Rotates the camera around its focus point by \p dx degrees on
/// the x-axis and \p dy degrees on the y-axis. If \p rotateDirection
/// is \c true then the direction will also be rotated so that the
/// camera remains pointed towards its focus point.
///
/// Equivalent to:
/// \code
/// orbit(focus(), dx, dy, rotateDirection);
/// \endcode
void GraphicsCamera::orbit(float dx, float dy, bool rotateDirection)
{
    orbit(d->focus, dx, dy, rotateDirection);
}

/// Rotates the camera around \p point by \p dx degrees on the x-axis
/// and \p dy degrees on the y-axis. If \p rotateDirection is \c true
/// the direction will also be rotated so that the camera remains
/// pointed toward \p point.
void GraphicsCamera::orbit(const Point3f &point, float dx, float dy, bool rotateDirection)
{
    moveBy(Vector3f(-point.x(), -point.y(), -point.z()));

    Vector3f up = upVector();
    Vector3f right = up.cross(direction()).normalized();

    rotate(up, dx, rotateDirection);
    rotate(right, dy, rotateDirection);

    moveBy(Vector3f(point.x(), point.y(), point.z()));
}

// --- Orientation --------------------------------------------------------- //
/// Sets the camera's direction.
void GraphicsCamera::setDirection(const Vector3f &direction)
{
    d->direction = direction.normalized();
    setChanged(true);
}

/// Returns the camera's direction.
Vector3f GraphicsCamera::direction() const
{
    return d->direction;
}

/// Sets the camera's focus to \p point.
void GraphicsCamera::setFocus(const Point3f &point)
{
    d->focus = point;
}

/// Returns the camera's focus.
Point3f GraphicsCamera::focus() const
{
    return d->focus;
}

/// Sets the camera's focus to \p point and rotates the camera's
/// direction to look at \p point.
///
/// Equivalent to:
/// \code
/// setFocus(point);
/// setDirection(point - position());
/// \endcode
void GraphicsCamera::lookAt(const Point3f &point)
{
    setFocus(point);
    setDirection(point - d->position);
}

/// Sets the camera's up vector.
void GraphicsCamera::setUpVector(const Vector3f &upVector)
{
    d->upVector = upVector.normalized();
    setChanged(true);
}

/// Returns the camera's up vector.
Vector3f GraphicsCamera::upVector() const
{
    return d->upVector;
}

/// Tilts the camera by \p angle degrees.
void GraphicsCamera::tilt(float angle)
{
    setUpVector(chemkit::geometry::rotate(upVector(), direction(), angle));
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsCamera::setChanged(bool changed)
{
    d->changed = changed;
}

bool GraphicsCamera::changed() const
{
    return d->changed;
}

} // end chemkit namespace
