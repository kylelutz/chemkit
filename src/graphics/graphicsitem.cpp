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

#include "graphicsitem.h"

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
