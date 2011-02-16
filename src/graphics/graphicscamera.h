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

#include "graphicspoint.h"
#include "graphicsvector.h"

namespace chemkit {

class GraphicsView;
class GraphicsCameraPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCamera
{
    public:
        // construction and destruction
        GraphicsCamera();
        GraphicsCamera(const GraphicsPoint &position);
        GraphicsCamera(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        ~GraphicsCamera();

        // properties
        GraphicsView* view() const;

        // position
        void setPosition(const GraphicsPoint &position);
        void setPosition(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        GraphicsPoint position() const;
        GraphicsFloat x() const;
        GraphicsFloat y() const;
        GraphicsFloat z() const;
        void moveTo(const GraphicsPoint &position);
        void moveTo(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        void moveBy(const GraphicsVector &vector);
        void moveBy(GraphicsFloat dx, GraphicsFloat dy, GraphicsFloat dz);
        void moveBy(GraphicsFloat distance, const GraphicsVector &direction);
        void moveFoward(GraphicsFloat distance);
        void moveBackward(GraphicsFloat distance);
        void rotate(const GraphicsVector &axis, GraphicsFloat angle, bool rotateDirection = true);
        void orbit(const GraphicsPoint &point, GraphicsFloat dx, GraphicsFloat dy, bool rotateDirection = true);

        // orientation
        void setDirection(const GraphicsVector &direction);
        GraphicsVector direction() const;
        void lookAt(const GraphicsPoint &position);
        void setUpVector(const GraphicsVector &upVector);
        GraphicsVector upVector() const;
        void tilt(GraphicsFloat angle);

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
