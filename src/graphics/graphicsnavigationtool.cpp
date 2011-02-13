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

#include "graphicsnavigationtool.h"

#include "graphicsview.h"

namespace chemkit {

// === GraphicsNavigationToolPrivate ======================================= //
class GraphicsNavigationToolPrivate
{
    public:
        bool mouseDown;
        QPoint lastPosition;
};

// === GraphicsNavigationTool ============================================== //
/// \class GraphicsNavigationTool graphicsnavigationtool.h chemkit/graphicsnavigationtool.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsNavigationTool class implements navigation
///        controls.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new navigation tool object.
GraphicsNavigationTool::GraphicsNavigationTool()
    : GraphicsTool(),
      d(new GraphicsNavigationToolPrivate)
{
    d->mouseDown = false;
}

/// Destroys the navigation tool object.
GraphicsNavigationTool::~GraphicsNavigationTool()
{
    delete d;
}

// --- Events -------------------------------------------------------------- //
void GraphicsNavigationTool::mousePressEvent(QMouseEvent *event)
{
    d->mouseDown = true;
    d->lastPosition = event->pos();
}

void GraphicsNavigationTool::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event);

    d->mouseDown = false;
}

void GraphicsNavigationTool::mouseMoveEvent(QMouseEvent *event)
{
    if(d->mouseDown){
        GraphicsCamera *camera = view()->camera();
        int dx = event->x() - d->lastPosition.x();
        int dy = event->y() - d->lastPosition.y();

        // left button
        if(event->buttons() & Qt::LeftButton){
            camera->orbit(GraphicsPoint(0, 0, 0), -dx, dy);
        }
        // right button
        if(event->buttons() & Qt::RightButton){
            camera->tilt(-dx);
            camera->moveFoward(dy);
        }

        d->lastPosition = event->pos();
        view()->update();
    }
}

void GraphicsNavigationTool::wheelEvent(QWheelEvent *event)
{
    GraphicsCamera *camera = view()->camera();

    if(event->delta() > 0){
        camera->moveFoward(5);
    }
    else{
        camera->moveBackward(5);
    }

    view()->update();
}

} // end chemkit namespace
