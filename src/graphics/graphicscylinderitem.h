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

#ifndef CHEMKIT_GRAPHICSCYLINDERITEM_H
#define CHEMKIT_GRAPHICSCYLINDERITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class GraphicsCylinderItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCylinderItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsCylinderItem(const GraphicsPoint &top, const GraphicsPoint &bottom, GraphicsFloat radius);
        ~GraphicsCylinderItem();

        // properties
        void setTop(const GraphicsPoint &top);
        GraphicsPoint top() const;
        void setBottom(const GraphicsPoint &bottom);
        GraphicsPoint bottom() const;
        void setRadius(GraphicsFloat radius);
        GraphicsFloat radius() const;
        void setColor(const QColor &color);
        QColor color() const;

        // drawing
        virtual void paint(GraphicsPainter *painter);

    private:
        GraphicsCylinderItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSCYLINDERITEM_H
