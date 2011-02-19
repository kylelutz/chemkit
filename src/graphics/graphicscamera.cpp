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

#include "graphicscamera.h"

#include "graphicsquaternion.h"

namespace chemkit {

// === GraphicsCameraPrivate =============================================== //
class GraphicsCameraPrivate
{
    public:
        GraphicsPoint position;
        GraphicsVector direction;
        GraphicsVector upVector;
        GraphicsPoint focus;
        GraphicsView *view;
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
    d->position = GraphicsPoint(0, 0, 0);
    d->direction = -GraphicsVector::Z();
    d->upVector = GraphicsVector::Y();
    d->focus = GraphicsPoint(0, 0, 0);
    d->changed = true;
    d->view = 0;
}

/// Creates a new graphics camera object at \p position.
GraphicsCamera::GraphicsCamera(const GraphicsPoint &position)
    : d(new GraphicsCameraPrivate)
{
    d->position = position;
    d->direction = -GraphicsVector::Z();
    d->upVector = GraphicsVector::Y();
    d->focus = GraphicsPoint(0, 0, 0);
    d->changed = true;
    d->view = 0;
}

/// \overload
GraphicsCamera::GraphicsCamera(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
    : d(new GraphicsCameraPrivate)
{
    d->position = GraphicsPoint(x, y, z);
    d->direction = -GraphicsVector::Z();
    d->upVector = GraphicsVector::Y();
    d->focus = GraphicsPoint(0, 0, 0);
    d->changed = true;
    d->view = 0;
}

/// Destroys the camera object.
GraphicsCamera::~GraphicsCamera()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the graphics view that the camera belongs to.
GraphicsView* GraphicsCamera::view() const
{
    return d->view;
}

// --- Position ------------------------------------------------------------ //
/// Sets the position of the camera to \p position.
void GraphicsCamera::setPosition(const GraphicsPoint &position)
{
    d->position = position;
    setChanged(true);
}

/// Sets the position of the camera to (\p x, \p y, \p z).
void GraphicsCamera::setPosition(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
{
    setPosition(GraphicsPoint(x, y, z));
}

/// Returns the position of the camera.
GraphicsPoint GraphicsCamera::position() const
{
    return d->position;
}

/// Returns the x component of the camera's position.
GraphicsFloat GraphicsCamera::x() const
{
    return position().x();
}

/// Returns the y component of the camera's position.
GraphicsFloat GraphicsCamera::y() const
{
    return position().y();
}

/// Returns the z component of the camera's position.
GraphicsFloat GraphicsCamera::z() const
{
    return position().z();
}

/// Sets the camera's position to \p position.
void GraphicsCamera::moveTo(const GraphicsPoint &position)
{
    setPosition(position);
}

/// Sets the camera's postion to (\p x, \p y, \p z).
void GraphicsCamera::moveTo(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z)
{
    setPosition(x, y, z);
}

/// Moves the camera by \p vector.
void GraphicsCamera::moveBy(const GraphicsVector &vector)
{
    setPosition(position().movedBy(vector));
}

/// Moves the camera by (\p dx, \p dy, \p dz).
void GraphicsCamera::moveBy(GraphicsFloat dx, GraphicsFloat dy, GraphicsFloat dz)
{
    moveBy(GraphicsVector(dx, dy, dz));
}

/// Moves the camera by \p distance in \p direction.
void GraphicsCamera::moveBy(GraphicsFloat distance, const GraphicsVector &direction)
{
    setPosition(position().movedBy(distance, direction));
}

/// Moves the camera fowards by \p distance.
void GraphicsCamera::moveFoward(GraphicsFloat distance)
{
    moveBy(distance, direction());
}

/// Moves the camera backwards by \p distance.
void GraphicsCamera::moveBackward(GraphicsFloat distance)
{
    moveBy(-distance, direction());
}

/// Rotates the camera around \p axis by \p angle degrees. If
/// \p rotateDirection is \c true the camera's direction will also
/// be rotated.
void GraphicsCamera::rotate(const GraphicsVector &axis, GraphicsFloat angle, bool rotateDirection)
{
    setPosition(GraphicsQuaternion::rotate(position(), axis, angle));

    if(rotateDirection){
        setDirection(GraphicsQuaternion::rotate(direction(), axis, angle));
        setUpVector(GraphicsQuaternion::rotate(upVector(), axis, angle));
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
void GraphicsCamera::orbit(GraphicsFloat dx, GraphicsFloat dy, bool rotateDirection)
{
    orbit(d->focus, dx, dy, rotateDirection);
}

/// Rotates the camera around \p point by \p dx degrees on the x-axis
/// and \p dy degrees on the y-axis. If \p rotateDirection is \c true
/// the direction will also be rotated so that the camera remains
/// pointed toward \p point.
void GraphicsCamera::orbit(const GraphicsPoint &point, GraphicsFloat dx, GraphicsFloat dy, bool rotateDirection)
{
    moveBy(GraphicsVector(-point.x(), -point.y(), -point.z()));

    GraphicsVector up = upVector();
    GraphicsVector right = up.cross(direction()).normalized();

    rotate(up, dx, rotateDirection);
    rotate(right, dy, rotateDirection);

    moveBy(GraphicsVector(point.x(), point.y(), point.z()));
}

// --- Orientation --------------------------------------------------------- //
/// Sets the camera's direction.
void GraphicsCamera::setDirection(const GraphicsVector &direction)
{
    d->direction = direction.normalized();
    setChanged(true);
}

/// Returns the camera's direction.
GraphicsVector GraphicsCamera::direction() const
{
    return d->direction;
}

/// Sets the camera's focus to \p point.
void GraphicsCamera::setFocus(const GraphicsPoint &point)
{
    d->focus = point;
}

/// Returns the camera's focus.
GraphicsPoint GraphicsCamera::focus() const
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
void GraphicsCamera::lookAt(const GraphicsPoint &point)
{
    setFocus(point);
    setDirection(point - d->position);
}

/// Sets the camera's up vector.
void GraphicsCamera::setUpVector(const GraphicsVector &upVector)
{
    d->upVector = upVector.normalized();
    setChanged(true);
}

/// Returns the camera's up vector.
GraphicsVector GraphicsCamera::upVector() const
{
    return d->upVector;
}

/// Tilts the camera by \p angle degrees.
void GraphicsCamera::tilt(GraphicsFloat angle)
{
    setUpVector(GraphicsQuaternion::rotate(upVector(), direction(), angle));
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsCamera::setView(GraphicsView *view)
{
    d->view = view;
}

void GraphicsCamera::setChanged(bool changed)
{
    d->changed = changed;
}

bool GraphicsCamera::changed() const
{
    return d->changed;
}

} // end chemkit namespace
