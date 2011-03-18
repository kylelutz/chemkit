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

#ifndef CHEMKIT_GRAPHICSPOINT_H
#define CHEMKIT_GRAPHICSPOINT_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/genericpoint.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsPoint : public GenericPoint<GraphicsFloat>
{
    public:
        // construction and destruction
        GraphicsPoint();
        GraphicsPoint(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        GraphicsPoint(const GenericPoint<float> &point);
        GraphicsPoint(const GenericPoint<double> &point);
        GraphicsPoint(const StaticVector<float, 3> &point);
        GraphicsPoint(const StaticVector<double, 3> &point);

        // properties
        Point3 toPoint3() const;
};

} // end chemkit namespace

#include "graphicspoint-inline.h"

#endif // CHEMKIT_GRAPHICSPOINT_H
