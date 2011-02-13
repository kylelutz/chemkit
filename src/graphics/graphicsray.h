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

#ifndef CHEMKIT_GRAPHICSRAY_H
#define CHEMKIT_GRAPHICSRAY_H

#include "graphics.h"

#include "graphicspoint.h"
#include "graphicsvector.h"

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsRay
{
    public:
        // construction and destruction
        GraphicsRay();
        GraphicsRay(const GraphicsPoint &origin, const GraphicsVector &direction);
        GraphicsRay(const GraphicsPoint &origin, const GraphicsPoint &point);
        ~GraphicsRay();

        // properties
        void setOrigin(const GraphicsPoint &origin);
        GraphicsPoint origin() const;
        void setDirection(const GraphicsVector &direction);
        GraphicsVector direction() const;

        // geometry
        GraphicsPoint pointAt(GraphicsFloat distance) const;
        bool intersectsSphere(const GraphicsPoint &center, GraphicsFloat radius, GraphicsFloat *distance = 0) const;
        bool intersectsCylinder(const GraphicsPoint &a, const GraphicsPoint &b, GraphicsFloat radius, GraphicsFloat *distance = 0) const;

    private:
        GraphicsPoint m_origin;
        GraphicsVector m_direction;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSRAY_H
