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

#ifndef CHEMKIT_GRAPHICSSCENE_H
#define CHEMKIT_GRAPHICSSCENE_H

#include "graphics.h"

#include "graphicsray.h"

namespace chemkit {

class GraphicsItem;
class GraphicsView;
class GraphicsScenePrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsScene
{
    public:
        // construction and destruction
        GraphicsScene();
        ~GraphicsScene();

        // properties
        int size() const;
        bool isEmpty() const;
        QList<GraphicsView *> views() const;

        // items
        void addItem(GraphicsItem *item);
        bool removeItem(GraphicsItem *item);
        bool deleteItem(GraphicsItem *item);
        GraphicsItem* item(int index) const;
        GraphicsItem* item(const GraphicsRay &ray) const;
        QList<GraphicsItem *> items() const;
        QList<GraphicsItem *> items(const GraphicsRay &ray, bool sorted = true) const;
        int itemCount() const;

    private:
        // internal methods
        void addView(GraphicsView *view);
        void removeView(GraphicsView *view);

        friend class GraphicsView;

    private:
        GraphicsScenePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSSCENE_H
