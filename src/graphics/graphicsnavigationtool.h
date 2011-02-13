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

#ifndef CHEMKIT_GRAPHICSNAVIGATIONTOOL_H
#define CHEMKIT_GRAPHICSNAVIGATIONTOOL_H

#include "graphics.h"

#include "graphicstool.h"

namespace chemkit {

class GraphicsNavigationToolPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsNavigationTool : public virtual GraphicsTool
{
    public:
        // construction and destruction
        GraphicsNavigationTool();
        ~GraphicsNavigationTool();

        // events
        virtual void mousePressEvent(QMouseEvent *event);
        virtual void mouseReleaseEvent(QMouseEvent *event);
        virtual void mouseMoveEvent(QMouseEvent *event);
        virtual void wheelEvent(QWheelEvent *event);

    private:
        GraphicsNavigationToolPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSNAVIGATIONTOOL_H
