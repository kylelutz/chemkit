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

#ifndef CHEMKIT_GRAPHICSOVERLAY_H
#define CHEMKIT_GRAPHICSOVERLAY_H

#include "graphics.h"

#include <chemkit/point3.h>

#include "graphicsitem.h"

namespace chemkit {

class GraphicsView;
class GraphicsOverlayPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsOverlay : public QGraphicsScene
{
    Q_OBJECT

    public:
        // construction and destruction
        GraphicsOverlay();
        ~GraphicsOverlay();

        // items
        void removeItem(QGraphicsItem *item);

        // binding
        void bindItemTo(QGraphicsItem *item, const Point3f &position);
        void bindItemTo(QGraphicsItem *item, GraphicsItem *sceneItem);
        void removeBinding(QGraphicsItem *item);
        void updateBindings(GraphicsView *view);

    private:
        GraphicsOverlayPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSOVERLAY_H
