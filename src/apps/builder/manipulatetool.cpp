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

#include "manipulatetool.h"

#include <set>

#include <chemkit/atom.h>
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
            editor()->setAtomPosition(atom, view()->unproject(event->x(), event->y(), atom->position().cast<float>()).cast<chemkit::Real>());
        }
        else if(event->buttons() & Qt::RightButton){
            int dy = event->y() - m_lastPosition.y();

            chemkit::Point3f position = atom->position().cast<float>();
            position += -view()->camera()->direction().normalized() * (dy * 0.1);
            editor()->setAtomPosition(atom, position.cast<chemkit::Real>());
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
    std::vector<chemkit::Atom *> newAtoms = editor()->paste();
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
    std::set<chemkit::Atom *> selection;

    for(int x = rect.x(); x < rect.right(); x += 5){
        for(int y = rect.y(); y < rect.bottom(); y += 5){
            chemkit::GraphicsItem *item = view()->itemAt(x, y);
            if(item && item->type() == chemkit::GraphicsItem::AtomItem){
                const chemkit::Atom *atom = static_cast<chemkit::GraphicsAtomItem *>(item)->atom();
                selection.insert(const_cast<chemkit::Atom *>(atom));
            }
        }
    }

    if(selection.empty()){
        clearSelection();
    }
    else{
        m_selection = std::vector<chemkit::Atom *>(selection.begin(), selection.end());
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
        chemkit::Vector3f delta = view()->unproject(finalPosition.x(), finalPosition.y(), atom->position().cast<float>()) -
                                  view()->unproject(initialPosition.x(), initialPosition.y(), atom->position().cast<float>());
        editor()->setAtomPosition(atom, atom->position() + delta.cast<chemkit::Real>());
    }

    view()->update();
}
