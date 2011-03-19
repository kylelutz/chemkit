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

#include "manipulatetool.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>
#include <chemkit/graphicscamera.h>
#include <chemkit/graphicsoverlay.h>

// --- Construction and Destruction ---------------------------------------- //
ManipulateTool::ManipulateTool(BuilderWindow *builder)
    : BuilderTool(builder)
{
    m_selectedItem = 0;
    m_selectionOverlayItem = 0;
    m_hasSelection = false;
}

ManipulateTool::~ManipulateTool()
{
}

// --- Events -------------------------------------------------------------- //
void ManipulateTool::mousePressEvent(QMouseEvent *event)
{
    m_initialPosition = event->posF();
    m_lastPosition = event->posF();

    if(m_hasSelection && view()->overlay()->itemAt(event->pos()) == m_selectionOverlayItem){
        setState(MovingSelection);
        builder()->beginMoleculeEdit();
    }
    else{
        clearSelection();

        chemkit::GraphicsItem *item = view()->itemAt(event->x(), event->y());

        if(item && item->type() == chemkit::GraphicsItem::AtomItem){
            m_selectedItem = item;
            builder()->beginMoleculeEdit();
            setState(MovingAtom);
        }
        else{
            m_selectedItem = 0;

            m_selectionOverlayItem = new QGraphicsRectItem(event->x(), event->y(), 0, 0);
            m_selectionOverlayItem->setPen(QPen(QBrush(Qt::white), 2, Qt::DashLine));
            m_selectionOverlayItem->setBrush(QBrush(Qt::green));
            m_selectionOverlayItem->setOpacity(0.55);
            view()->overlay()->addItem(m_selectionOverlayItem);
            setState(Selecting);
        }
    }
}

void ManipulateTool::mouseMoveEvent(QMouseEvent *event)
{
    if(state() == MovingSelection){
        QPointF delta = event->pos() - m_lastPosition;
        moveSelectionBy(delta.x(), delta.y());
    }
    else if(state() == Selecting){
        m_selectionOverlayItem->setRect(QRectF(m_initialPosition, event->pos()).normalized());
    }
    else if(state() == MovingAtom){
        chemkit::GraphicsAtomItem *atomItem = static_cast<chemkit::GraphicsAtomItem *>(m_selectedItem);
        chemkit::Atom *atom = const_cast<chemkit::Atom *>(atomItem->atom());

        if(event->buttons() & Qt::LeftButton){
            editor()->setAtomPosition(atom, view()->unproject(event->x(), event->y(), atom->position()));
        }
        else if(event->buttons() & Qt::RightButton){
            int dy = event->y() - m_lastPosition.y();

            chemkit::Point3f position = atom->position();
            position.moveBy(dy * 0.1, -view()->camera()->direction());
            editor()->setAtomPosition(atom, position);
        }
    }

    m_lastPosition = event->posF();
    view()->update();
}

void ManipulateTool::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event);

    if(state() == MovingAtom){
        builder()->endMoleculeEdit();
    }
    else if(state() == MovingSelection){
        builder()->endMoleculeEdit();
    }
    else if(state() == Selecting){
        setSelection(m_selectionOverlayItem->rect().toRect().normalized());
    }

    m_selectedItem = 0;
    view()->update();
}

void ManipulateTool::toolChanged(const GraphicsTool *tool)
{
    Q_UNUSED(tool);

    clearSelection();
}

void ManipulateTool::cut()
{
    editor()->cut(m_selection);
}

void ManipulateTool::copy()
{
    editor()->copy(m_selection);
}

void ManipulateTool::paste()
{
    builder()->beginMoleculeEdit();
    QList<chemkit::Atom *> newAtoms = editor()->paste();
    m_selection = newAtoms;
    moveSelectionBy(30, -30);
    builder()->endMoleculeEdit();
}

void ManipulateTool::del()
{
    builder()->beginMoleculeEdit();

    foreach(chemkit::Atom *atom, m_selection){
        editor()->removeAtom(atom);
    }

    builder()->endMoleculeEdit();

    clearSelection();
}

// --- Internal Methods ---------------------------------------------------- //
void ManipulateTool::setState(State state)
{
    if(m_state == state)
        return;

    m_state = state;
}

ManipulateTool::State ManipulateTool::state() const
{
    return m_state;
}

void ManipulateTool::setSelection(const QRect &rect)
{
    QSet<chemkit::Atom *> selection;

    for(int x = rect.x(); x < rect.right(); x += 5){
        for(int y = rect.y(); y < rect.bottom(); y += 5){
            chemkit::GraphicsItem *item = view()->itemAt(x, y);
            if(item && item->type() == chemkit::GraphicsItem::AtomItem){
                const chemkit::Atom *atom = static_cast<chemkit::GraphicsAtomItem *>(item)->atom();
                selection.insert(const_cast<chemkit::Atom *>(atom));
            }
        }
    }

    if(selection.isEmpty()){
        clearSelection();
    }
    else{
        m_selection = selection.toList();
        m_hasSelection = true;

        setCanCut(true);
        setCanCopy(true);
        setCanDelete(true);
    }
}

void ManipulateTool::clearSelection()
{
    m_selection.clear();

    if(m_selectionOverlayItem){
        view()->overlay()->removeItem(m_selectionOverlayItem);
        delete m_selectionOverlayItem;
        m_selectionOverlayItem = 0;
        view()->update();
    }

    m_hasSelection = false;

    setCanCut(false);
    setCanCopy(false);
    setCanDelete(false);

    editor()->clearCopyBuffer();
}

void ManipulateTool::moveSelectionBy(int x, int y)
{
    if(m_selectionOverlayItem)
        m_selectionOverlayItem->moveBy(x, y);

    QPointF initialPosition = m_lastPosition;
    QPointF finalPosition = m_lastPosition + QPointF(x, y);

    foreach(chemkit::Atom *atom, m_selection){
        chemkit::Vector3f delta = view()->unproject(finalPosition.x(), finalPosition.y(), atom->position()) -
                                  view()->unproject(initialPosition.x(), initialPosition.y(), atom->position());
        editor()->setAtomPosition(atom, atom->position() + delta);
    }

    view()->update();
}
