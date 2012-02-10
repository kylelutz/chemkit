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

#include "graphicsitem.h"

#include <chemkit/foreach.h>

#include "graphicsview.h"
#include "graphicsscene.h"
#include "graphicsmaterial.h"
#include "graphicstransform.h"

namespace chemkit {

// === GraphicsItemPrivate ================================================= //
class GraphicsItemPrivate
{
public:
    int type;
    bool visible;
    GraphicsScene *scene;
    float opacity;
    GraphicsMaterial *material;
    GraphicsTransform transform;
};

// === GraphicsItem ======================================================== //
/// \class GraphicsItem graphicsitem.h chemkit/graphicsitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsItem class is the base class for all graphics
///        items.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics item object.
GraphicsItem::GraphicsItem(int type)
    : d(new GraphicsItemPrivate)
{
    d->type = type;
    d->visible = true;
    d->scene = 0;
    d->opacity = 1.0;
    d->material = new GraphicsMaterial;
    d->transform = GraphicsTransform::identity();
}

/// Destroys the graphics item object.
GraphicsItem::~GraphicsItem()
{
    if(d->scene){
        d->scene->removeItem(this);
    }

    delete d->material;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
int GraphicsItem::type() const
{
    return d->type;
}

/// Sets the item's visibility.
void GraphicsItem::setVisible(bool visible)
{
    d->visible = visible;
    itemChanged(ItemVisiblityChanged);
}

/// Returns \c true if the item is visible.
bool GraphicsItem::isVisible() const
{
    return d->visible;
}

/// Sets the item's visibility to \c true.
void GraphicsItem::show()
{
    setVisible(true);
}

/// Sets the item's visibility to \c false.
void GraphicsItem::hide()
{
    setVisible(false);
}

/// Returns the \p scene the item is contained in or \c 0 if the
/// item is not in a scene.
GraphicsScene* GraphicsItem::scene() const
{
    return d->scene;
}

void GraphicsItem::setOpacity(float opacity)
{
    d->opacity = opacity;

    itemChanged(ItemOpacityChanged);
}

float GraphicsItem::opacity() const
{
    return d->opacity;
}

bool GraphicsItem::isOpaque() const
{
    return opacity() > 0.99;
}

bool GraphicsItem::isTransparent() const
{
    return opacity() < 0.01;
}

bool GraphicsItem::isTranslucent() const
{
    return !(isOpaque() || isTransparent());
}

// --- Material ------------------------------------------------------------ //
void GraphicsItem::setMaterial(GraphicsMaterial *material)
{
    delete d->material;

    d->material = material;
}

GraphicsMaterial* GraphicsItem::material() const
{
    return d->material;
}

// --- Geometry ------------------------------------------------------------ //
void GraphicsItem::setTransform(const GraphicsTransform &transform)
{
    d->transform = transform;
}

GraphicsTransform GraphicsItem::transform() const
{
    return d->transform;
}

void GraphicsItem::translate(const Vector3f &vector)
{
    d->transform *= GraphicsTransform::translation(vector);
}

void GraphicsItem::translate(float x, float y, float z)
{
    translate(Vector3f(x, y, z));
}

void GraphicsItem::rotate(const Vector3f &axis, const float angle)
{
    d->transform *= GraphicsTransform::rotation(axis, angle);
}

/// Returns the axis-aligned bounding box for the item.
GraphicsBoundingBox GraphicsItem::boundingBox() const
{
    return GraphicsBoundingBox();
}

bool GraphicsItem::intersects(const GraphicsRay &ray, float *distance) const
{
    Q_UNUSED(ray);
    Q_UNUSED(distance);

    return false;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsItem::paint(GraphicsPainter *painter)
{
    Q_UNUSED(painter);
}

void GraphicsItem::update()
{
    if(scene()){
        foreach(GraphicsView *view, scene()->views()){
            view->update();
        }
    }
}

// --- Events -------------------------------------------------------------- //
void GraphicsItem::itemChanged(ItemChange change)
{
    Q_UNUSED(change);
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsItem::setScene(GraphicsScene *scene)
{
    d->scene = scene;
    itemChanged(ItemSceneChanged);
}

} // end chemkit namespace
