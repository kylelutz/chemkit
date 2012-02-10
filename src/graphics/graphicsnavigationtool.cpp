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

#include "graphicsnavigationtool.h"

#include "graphicsview.h"
#include "graphicscamera.h"

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
        const boost::shared_ptr<GraphicsCamera> &camera = view()->camera();
        int dx = event->x() - d->lastPosition.x();
        int dy = event->y() - d->lastPosition.y();

        // left button
        if(event->buttons() & Qt::LeftButton){
            camera->orbit(-dx, dy);
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
    const boost::shared_ptr<GraphicsCamera> &camera = view()->camera();

    if(event->delta() > 0){
        camera->moveFoward(5);
    }
    else{
        camera->moveBackward(5);
    }

    view()->update();
}

} // end chemkit namespace
