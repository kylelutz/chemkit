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

#include "graphicsoverlay.h"

#include "graphicsview.h"

namespace chemkit {

// === GraphicsOverlayPrivate ============================================== //
class GraphicsOverlayPrivate
{
    public:
        QHash<QGraphicsItem *, Point3f> bindings;
};

// === GraphicsOverlay ===================================================== //
/// \class GraphicsOverlay graphicsoverlay.h chemkit/graphicsoverlay.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsOverlay class represents an overlay on a
///        graphics view.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new overlay.
GraphicsOverlay::GraphicsOverlay()
    : QGraphicsScene(),
      d(new GraphicsOverlayPrivate)
{
}

/// Destroys the overlay object.
GraphicsOverlay::~GraphicsOverlay()
{
    delete d;
}

// --- Items --------------------------------------------------------------- //
void GraphicsOverlay::removeItem(QGraphicsItem *item)
{
    QGraphicsScene::removeItem(item);

    if(d->bindings.contains(item)){
        d->bindings.remove(item);
    }
}

// --- Binding ------------------------------------------------------------- //
void GraphicsOverlay::bindItemTo(QGraphicsItem *item, const Point3f &position)
{
    d->bindings[item] = position;
}

void GraphicsOverlay::removeBinding(QGraphicsItem *item)
{
    d->bindings.remove(item);
}

void GraphicsOverlay::updateBindings(GraphicsView *view)
{
    foreach(QGraphicsItem *item, d->bindings.keys()){
        QPointF itemPos = view->project(d->bindings[item]);

        // center item
        QRectF rect = item->boundingRect();
        itemPos.setX(itemPos.x() - rect.width()/2);
        itemPos.setY(itemPos.y() - rect.height()/2);
        item->setPos(itemPos);
    }
}

} // end chemkit namespace
