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

#ifndef CHEMKIT_GRAPHICSQUATERNION_H
#define CHEMKIT_GRAPHICSQUATERNION_H

#include "graphics.h"

#include <chemkit/genericquaternion.h>

#include "graphicspoint.h"
#include "graphicsvector.h"

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsQuaternion : public GenericQuaternion<GraphicsFloat>
{
    public:
        // construction and destruction
        GraphicsQuaternion(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z, GraphicsFloat r);
        GraphicsQuaternion(const GraphicsPoint &point, GraphicsFloat r);
        GraphicsQuaternion(const GraphicsVector &vector, GraphicsFloat r);
        GraphicsQuaternion(const GenericQuaternion<GraphicsFloat> &quaternion);
        GraphicsQuaternion(const StaticVector<GraphicsFloat, 4> &quaternion);

        // properties
        GraphicsPoint toPoint() const;
        GraphicsVector toVector() const;

        // static methods
        static GraphicsQuaternion rotation(const GraphicsVector &axis, GraphicsFloat angle);
        static GraphicsQuaternion rotationRadians(const GraphicsVector &axis, GraphicsFloat angle);
        static GraphicsPoint rotate(const GraphicsPoint &point, const GraphicsVector &axis, GraphicsFloat angle);
        static GraphicsPoint rotateRadians(const GraphicsPoint &point, const GraphicsVector &axis, GraphicsFloat angle);
        static GraphicsVector rotate(const GraphicsVector &vector, const GraphicsVector &axis, GraphicsFloat angle);
        static GraphicsVector rotateRadians(const GraphicsVector &vector, const GraphicsVector &axis, GraphicsFloat angle);
};

} // end chemkit namespace

#include "graphicsquaternion-inline.h"

#endif // CHEMKIT_GRAPHICSQUATERNION_H
