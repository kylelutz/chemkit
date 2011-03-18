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

#ifndef CHEMKIT_GRAPHICSRINGITEM_H
#define CHEMKIT_GRAPHICSRINGITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class Ring;
class GraphicsRingItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsRingItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsRingItem(const Ring *ring = 0);
        ~GraphicsRingItem();

        // properties
        void setRing(const Ring *ring);
        const Ring* ring() const;
        void setColor(const QColor &color);
        QColor color() const;

        // drawing
        virtual void paint(GraphicsPainter *painter);

    private:
        Point3g ringCenter(const Ring *ring) const;

    private:
        GraphicsRingItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSRINGITEM_H
