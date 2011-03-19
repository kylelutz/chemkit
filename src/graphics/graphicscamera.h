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

#include "point3g.h"
#include "vector3g.h"

namespace chemkit {

class GraphicsView;
class GraphicsCameraPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCamera
{
    public:
        // construction and destruction
        GraphicsCamera();
        GraphicsCamera(const Point3g &position);
        GraphicsCamera(float x, float y, float z);
        ~GraphicsCamera();

        // properties
        GraphicsView* view() const;

        // position
        void setPosition(const Point3g &position);
        void setPosition(float x, float y, float z);
        Point3g position() const;
        float x() const;
        float y() const;
        float z() const;
        void moveTo(const Point3g &position);
        void moveTo(float x, float y, float z);
        void moveBy(const Vector3g &vector);
        void moveBy(float dx, float dy, float dz);
        void moveBy(float distance, const Vector3g &direction);
        void moveFoward(float distance);
        void moveBackward(float distance);
        void rotate(const Vector3g &axis, float angle, bool rotateDirection = true);
        void orbit(float dx, float dy, bool rotateDirection = true);
        void orbit(const Point3g &point, float dx, float dy, bool rotateDirection = true);

        // orientation
        void setDirection(const Vector3g &direction);
        Vector3g direction() const;
        void setFocus(const Point3g &point);
        Point3g focus() const;
        void lookAt(const Point3g &point);
        void setUpVector(const Vector3g &upVector);
        Vector3g upVector() const;
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
