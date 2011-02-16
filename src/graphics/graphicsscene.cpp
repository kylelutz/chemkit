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

#include "graphicsscene.h"

#include "graphicsitem.h"
#include "graphicsview.h"

namespace chemkit {

// === GraphicsScenePrivate ================================================ //
class GraphicsScenePrivate
{
    public:
        QList<GraphicsItem *> items;
        QList<GraphicsView *> views;
};

// === GraphicsScene ======================================================= //
/// \class GraphicsScene graphicsscene.h chemkit/graphicsscene.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsScene class contains graphics items.
///
/// The GraphicsScene class contains and organizes GraphicsItems.
///
/// To display a graphics scene use the GraphicsView class.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics scene.
GraphicsScene::GraphicsScene()
    : d(new GraphicsScenePrivate)
{
}

/// Destroys the graphics scene.
GraphicsScene::~GraphicsScene()
{
    foreach(GraphicsItem *item, d->items){
        deleteItem(item);
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of items in the scene.
int GraphicsScene::size() const
{
    return itemCount();
}

/// Returns \c true if the scene contains no items.
bool GraphicsScene::isEmpty() const
{
    return size() == 0;
}

/// Returns a list of views that show the scene.
QList<GraphicsView *> GraphicsScene::views() const
{
    return d->views;
}

// --- Items --------------------------------------------------------------- //
/// Adds \p item to the scene.
///
/// The scene takes ownership of the item.
void GraphicsScene::addItem(GraphicsItem *item)
{
    d->items.append(item);
    item->setScene(this);
}

/// Removes \p item from the scene. Returns \c true if the item was
/// found and removed successfully.
///
/// The ownership of item is passed to the caller.
bool GraphicsScene::removeItem(GraphicsItem *item)
{
    bool found = d->items.removeOne(item);

    if(found){
        item->setScene(0);
    }

    return found;
}

/// Removes \p item from the scene and deletes it. Returns \c true
/// if the item was found and removed successfully.
bool GraphicsScene::deleteItem(GraphicsItem *item)
{
    bool found = removeItem(item);

    if(found){
        delete item;
    }

    return found;
}

/// Returns the item at \p index.
GraphicsItem* GraphicsScene::item(int index) const
{
    return d->items.value(index, 0);
}

/// Returns the item that intersects \p ray.
GraphicsItem* GraphicsScene::item(const GraphicsRay &ray) const
{
    GraphicsItem *closestItem = 0;
    GraphicsFloat closestDistance = qInf();

    foreach(GraphicsItem *item, d->items){
        GraphicsFloat distance;

        if(item->intersects(ray, &distance)){
            if(!closestItem || distance < closestDistance){
                closestItem = item;
                closestDistance = distance;
            }
        }
    }

    return closestItem;
}

/// Returns a list of items in the scene.
QList<GraphicsItem *> GraphicsScene::items() const
{
    return d->items;
}

/// Returns a list of all items that intersect \p ray.
QList<GraphicsItem *> GraphicsScene::items(const GraphicsRay &ray, bool sorted) const
{
    QList<GraphicsItem *> items;

    if(sorted){
        QList<GraphicsFloat> distances;

        foreach(GraphicsItem *item, d->items){
            GraphicsFloat distance;

            if(item->intersects(ray, &distance)){
                int index;
                for(index = 0; index < items.size(); index++)
                    if(distance < distances[index])
                        break;

                items.insert(index, item);
                distances.insert(index, distance);
            }
        }
    }
    else{
        foreach(GraphicsItem *item, d->items){
            if(item->intersects(ray)){
                items.append(item);
            }
        }
    }

    return items;
}

/// Returns the number of items in the scene.
int GraphicsScene::itemCount() const
{
    return d->items.size();
}

// --- Internal Methods ---------------------------------------------------- //
void GraphicsScene::addView(GraphicsView *view)
{
    d->views.append(view);
}

void GraphicsScene::removeView(GraphicsView *view)
{
    d->views.removeOne(view);
}

} // end chemkit namespace
