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

#include "graphicssphereitem.h"

#include "graphicspainter.h"

namespace chemkit {

// === GraphicsSphereItemPrivate =========================================== //
class GraphicsSphereItemPrivate
{
    public:
        Point3g position;
        float radius;
        QColor color;
};

// === GraphicsSphereItem ================================================== //
/// \class GraphicsSphereItem graphicssphereitem.h chemkit/graphicssphereitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsSphereItem class represents a sphere.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new sphere item at \p position with \p radius.
GraphicsSphereItem::GraphicsSphereItem(const Point3g &position, float radius)
    : GraphicsItem(),
      d(new GraphicsSphereItemPrivate)
{
    setPosition(position);
    setRadius(radius);
    d->color = Qt::red;
}

/// Creates a new sphere item at \p (\p x, \p y, \p z) with \p radius.
GraphicsSphereItem::GraphicsSphereItem(float x, float y, float z, float radius)
    : GraphicsItem(),
      d(new GraphicsSphereItemPrivate)
{
    setPosition(x, y, z);
    setRadius(radius);
    d->color = Qt::red;
}

/// Destroys the sphere item.
GraphicsSphereItem::~GraphicsSphereItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the sphere's position to \p position.
void GraphicsSphereItem::setPosition(const Point3g &position)
{
    d->position = position;

    setTransform(GraphicsTransform::translation(position));
}

/// Sets the sphere's position to (\p x, \p y, \p z).
void GraphicsSphereItem::setPosition(float x, float y, float z)
{
    setPosition(Point3g(x, y, z));
}

/// Returns the sphere's position.
Point3g GraphicsSphereItem::position() const
{
    return d->position;
}

/// Sets the sphere's radius to \p radius.
void GraphicsSphereItem::setRadius(float radius)
{
    d->radius = radius;
}

/// Returns the sphere's radius.
float GraphicsSphereItem::radius() const
{
    return d->radius;
}

/// Sets the sphere's color to \p color.
void GraphicsSphereItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the sphere's color.
QColor GraphicsSphereItem::color() const
{
    return d->color;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsSphereItem::paint(GraphicsPainter *painter)
{
    painter->setColor(d->color);
    painter->drawSphere(d->radius);
}

} // end chemkit namespace
