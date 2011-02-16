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
void GraphicsTool::toolChanged(const GraphicsTool *tool)
{
    Q_UNUSED(tool);
}

} // end chemkit namespace
