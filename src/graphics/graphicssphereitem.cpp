/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "graphicssphereitem.h"

#include "graphicspainter.h"

namespace chemkit {

// === GraphicsSphereItemPrivate =========================================== //
class GraphicsSphereItemPrivate
{
public:
    Point3f position;
    float radius;
    QColor color;
};

// === GraphicsSphereItem ================================================== //
/// \class GraphicsSphereItem graphicssphereitem.h chemkit/graphicssphereitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsSphereItem class represents a sphere.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new sphere item at \p position with \p radius.
GraphicsSphereItem::GraphicsSphereItem(const Point3f &position, float radius)
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
void GraphicsSphereItem::setPosition(const Point3f &position)
{
    d->position = position;

    setTransform(GraphicsTransform::translation(position));
}

/// Sets the sphere's position to (\p x, \p y, \p z).
void GraphicsSphereItem::setPosition(float x, float y, float z)
{
    setPosition(Point3f(x, y, z));
}

/// Returns the sphere's position.
Point3f GraphicsSphereItem::position() const
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
