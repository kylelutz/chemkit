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

#ifndef CHEMKIT_GRAPHICSCAMERA_H
#define CHEMKIT_GRAPHICSCAMERA_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>

namespace chemkit {

class GraphicsView;
class GraphicsCameraPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCamera
{
    public:
        // construction and destruction
        GraphicsCamera();
        GraphicsCamera(const Point3f &position);
        GraphicsCamera(float x, float y, float z);
        ~GraphicsCamera();

        // properties
        GraphicsView* view() const;

        // position
        void setPosition(const Point3f &position);
        void setPosition(float x, float y, float z);
        Point3f position() const;
        float x() const;
        float y() const;
        float z() const;
        void moveTo(const Point3f &position);
        void moveTo(float x, float y, float z);
        void moveBy(const Vector3f &vector);
        void moveBy(float dx, float dy, float dz);
        void moveBy(float distance, const Vector3f &direction);
        void moveFoward(float distance);
        void moveBackward(float distance);
        void rotate(const Vector3f &axis, float angle, bool rotateDirection = true);
        void orbit(float dx, float dy, bool rotateDirection = true);
        void orbit(const Point3f &point, float dx, float dy, bool rotateDirection = true);

        // orientation
        void setDirection(const Vector3f &direction);
        Vector3f direction() const;
        void setFocus(const Point3f &point);
        Point3f focus() const;
        void lookAt(const Point3f &point);
        void setUpVector(const Vector3f &upVector);
        Vector3f upVector() const;
        void tilt(float angle);

    private:
        // internal methods
        void setView(GraphicsView *view);
        void setChanged(bool changed);
        bool changed() const;

        friend class GraphicsView;

    private:
        GraphicsCameraPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSCAMERA_H
