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

#include <chemkit/vector3.h>

#include "graphicsitem.h"

namespace chemkit {

class Bond;
class GraphicsBondItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsBondItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsBondItem(const Bond *bond, float radius = 0.15);
        ~GraphicsBondItem();

        // properties
        void setBond(const Bond *bond);
        const Bond* bond() const;
        void setRadius(float radius);
        float radius() const;
        void setMaximumRadius(float radius);
        float maximumRadius() const;
        void setNormal(const Vector3f &normal);
        Vector3f normal() const;
        void setAtomColored(bool atomColored);
        bool atomColored() const;
        void setBondOrderVisible(bool showBondOrder);
        bool bondOrderVisible() const;
        void setColor(const QColor &color);
        QColor color() const;
        void setAtomColors(const QColor &a, const QColor &b);
        QPair<QColor, QColor> atomColors();

        // intersection
        virtual bool intersects(const GraphicsRay &ray, float *distance = 0) const;

        // drawing
        void paint(GraphicsPainter *painter);

    private:
        GraphicsBondItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSBONDITEM_H
