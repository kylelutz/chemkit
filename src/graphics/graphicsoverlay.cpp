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

#include "graphicsoverlay.h"

#include <chemkit/foreach.h>

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
