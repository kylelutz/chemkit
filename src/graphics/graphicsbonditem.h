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

#ifndef CHEMKIT_GRAPHICSBONDITEM_H
#define CHEMKIT_GRAPHICSBONDITEM_H

#include "graphics.h"

#include "graphicsitem.h"
#include "graphicsvector.h"

namespace chemkit {

class Bond;
class GraphicsBondItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsBondItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsBondItem(const Bond *bond, GraphicsFloat radius = 0.15);
        ~GraphicsBondItem();

        // properties
        void setBond(const Bond *bond);
        const Bond* bond() const;
        void setRadius(GraphicsFloat radius);
        GraphicsFloat radius() const;
        void setMaximumRadius(GraphicsFloat radius);
        GraphicsFloat maximumRadius() const;
        void setNormal(const GraphicsVector &normal);
        GraphicsVector normal() const;
        void setAtomColored(bool atomColored);
        bool atomColored() const;
        void setBondOrderVisible(bool showBondOrder);
        bool bondOrderVisible() const;

        // intersection
        virtual bool intersects(const GraphicsRay &ray, GraphicsFloat *distance = 0) const;

        // drawing
        void paint(GraphicsPainter *painter);

    private:
        GraphicsBondItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSBONDITEM_H
