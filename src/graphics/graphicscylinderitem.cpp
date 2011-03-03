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

#include "graphicscylinderitem.h"

#include "graphicspainter.h"

namespace chemkit {

// === GraphicsCylinderItemPrivate ========================================= //
class GraphicsCylinderItemPrivate
{
    public:
        GraphicsPoint top;
        GraphicsPoint bottom;
        GraphicsFloat radius;
        QColor color;
};

// === GraphicsCylinderItem ================================================ //
/// \class GraphicsCylinderItem graphicscylinderitem.h chemkit/graphicscylinderitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsCylinderItem class represents a cylinder.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new cylinder item with positions \p top and \p bottom
/// and with \p radius.
GraphicsCylinderItem::GraphicsCylinderItem(const GraphicsPoint &top, const GraphicsPoint &bottom, GraphicsFloat radius)
    : GraphicsItem(),
      d(new GraphicsCylinderItemPrivate)
{
    d->top = top;
    d->bottom = bottom;
    d->radius = radius;
    d->color = Qt::red;
}

/// Destroys the cylinder item.
GraphicsCylinderItem::~GraphicsCylinderItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the top position of the cylinder to \p top.
void GraphicsCylinderItem::setTop(const GraphicsPoint &top)
{
    d->top = top;
}

/// Returns the top position of the cylinder.
GraphicsPoint GraphicsCylinderItem::top() const
{
    return d->top;
}

/// Sets the bottom position of the cylinder to \p bottom.
void GraphicsCylinderItem::setBottom(const GraphicsPoint &bottom)
{
    d->bottom = bottom;
}

/// Returns the bottom position of the cylinder.
GraphicsPoint GraphicsCylinderItem::bottom() const
{
    return d->bottom;
}

/// Sets the cylinder's radius to \p radius.
void GraphicsCylinderItem::setRadius(GraphicsFloat radius)
{
    d->radius = radius;
}

/// Returns the cylinder's radius.
GraphicsFloat GraphicsCylinderItem::radius() const
{
    return d->radius;
}

/// Sets the cylinder's color to \p color.
void GraphicsCylinderItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the cylinder's color.
QColor GraphicsCylinderItem::color() const
{
    return d->color;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsCylinderItem::paint(GraphicsPainter *painter)
{
    painter->setColor(d->color);

    // draw cylinder
    painter->drawCylinder(d->top, d->bottom, d->radius);

    // draw caps
    painter->drawCircle(d->top, d->radius, d->top - d->bottom);
    painter->drawCircle(d->bottom, d->radius, d->bottom - d->top);
}

} // end chemkit namespace
