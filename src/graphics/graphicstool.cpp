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

#include "graphicstool.h"

#include "graphicsview.h"

namespace chemkit {

// === GraphicsToolPrivate ================================================= //
class GraphicsToolPrivate
{
public:
    GraphicsView *view;
};

// === GraphicsTool ======================================================== //
/// \class GraphicsTool graphicstool.h chemkit/graphicstool.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsTool class handles graphics events for a
///        GraphicsView.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics tool.
GraphicsTool::GraphicsTool()
    : d(new GraphicsToolPrivate)
{
    d->view = 0;
}

/// Destroys the graphics tool.
GraphicsTool::~GraphicsTool()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the view to \p view.
void GraphicsTool::setView(GraphicsView *view)
{
    d->view = view;
}

/// Returns the view the tool is a part of. Returns \c 0 if the tool
/// is not in any view.
GraphicsView* GraphicsTool::view() const
{
    return d->view;
}

// --- Event Handling------------------------------------------------------- //
/// Handle a mouse press event.
void GraphicsTool::mousePressEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
}

/// Handle a mouse release event.
void GraphicsTool::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
}

/// Handle a mouse double click event.
void GraphicsTool::mouseDoubleClickEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
}

/// Handle a mouse move event.
void GraphicsTool::mouseMoveEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
}

/// Handle a mouse wheel event.
void GraphicsTool::wheelEvent(QWheelEvent *event)
{
    Q_UNUSED(event);
}

/// This method is called when the current tool in the view changes.
///
/// \see GraphicsView::setTool()
void GraphicsTool::toolChanged(const boost::shared_ptr<GraphicsTool> &tool)
{
    Q_UNUSED(tool);
}

} // end chemkit namespace
