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

#ifndef CHEMKIT_GRAPHICSPROTEINITEM_H
#define CHEMKIT_GRAPHICSPROTEINITEM_H

#include "graphics.h"

#include "graphicsitem.h"
#include "graphicsproteincoilitem.h"
#include "graphicsproteinhelixitem.h"
#include "graphicsproteinsheetitem.h"

namespace chemkit {

class Polymer;
class GraphicsProteinItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsProteinItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsProteinItem(const Polymer *polymer = 0);
        ~GraphicsProteinItem();

        // properties
        void setPolymer(const Polymer *polymer);
        const Polymer* polymer() const;
        void setSecondaryStructureVisible(bool visible);
        bool secondaryStructureVisible() const;
        void setCoilRadius(GraphicsFloat radius);
        GraphicsFloat coilRadius() const;
        void setHelixDisplayType(GraphicsProteinHelixItem::DisplayType type);
        GraphicsProteinHelixItem::DisplayType helixDisplayType() const;

        // painting
        virtual void paint(GraphicsPainter *painter);

    private:
        GraphicsProteinItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPROTEINITEM_H
