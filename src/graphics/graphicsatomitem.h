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

#ifndef CHEMKIT_GRAPHICSATOMITEM_H
#define CHEMKIT_GRAPHICSATOMITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class Atom;
class GraphicsAtomItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsAtomItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsAtomItem(const Atom *atom, float radius = 0.5);
        ~GraphicsAtomItem();

        // properties
        void setAtom(const Atom *atom);
        const Atom* atom() const;
        void setRadius(float radius);
        float radius() const;
        void setColor(const QColor &color);
        QColor color() const;

        // geometry
        bool intersects(const GraphicsRay &ray, float *distance = 0) const;

        // drawing
        virtual void paint(GraphicsPainter *painter);

    private:
        GraphicsAtomItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSATOMITEM_H
