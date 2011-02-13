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

#ifndef MANIPULATETOOL_H
#define MANIPULATETOOL_H

#include "buildertool.h"

class ManipulateTool : public BuilderTool
{
    public:
        // construction and destruction
        ManipulateTool(BuilderWindow *builder);
        ~ManipulateTool();

        // events
        virtual void mousePressEvent(QMouseEvent *event);
        virtual void mouseReleaseEvent(QMouseEvent *event);
        virtual void mouseMoveEvent(QMouseEvent *event);
        virtual void toolChanged(const GraphicsTool *tool);
        virtual void cut();
        virtual void copy();
        virtual void paste();
        virtual void del();

    private:
        enum State{
            MovingAtom,
            Selecting,
            MovingSelection
        };

        void setState(State state);
        State state() const;
        void setSelection(const QRect &rect);
        void clearSelection();
        void moveSelectionBy(int x, int y);

    private:
        chemkit::GraphicsItem *m_selectedItem;
        QGraphicsRectItem *m_selectionOverlayItem;
        QPointF m_initialPosition;
        QPointF m_lastPosition;
        QList<chemkit::Atom *> m_selection;
        bool m_hasSelection;
        enum State m_state;
};

#endif // MANIPULATETOOL_H
