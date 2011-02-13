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

#ifndef CHEMKIT_GRAPHICSCYLINDER_H
#define CHEMKIT_GRAPHICSCYLINDER_H

#include "graphics.h"

namespace chemkit {

class GraphicsRay;
class GraphicsVertexBuffer;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCylinder
{
    public:
        // construction and destruction
        GraphicsCylinder();
        GraphicsCylinder(GraphicsFloat radius, GraphicsFloat length);
        ~GraphicsCylinder();

        // properties
        void setRadius(GraphicsFloat radius);
        GraphicsFloat radius() const;
        void setLength(GraphicsFloat length);
        GraphicsFloat length() const;

        // intersection
        bool intersects(const GraphicsRay &ray, GraphicsFloat *distance = 0) const;

        // tesselation
        GraphicsVertexBuffer* tesselate(int slices, int stacks) const;

    private:
        GraphicsFloat m_radius;
        GraphicsFloat m_length;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSCYLINDER_H
