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

#ifndef CHEMKIT_GRAPHICSPROTEINHELIXITEM_H
#define CHEMKIT_GRAPHICSPROTEINHELIXITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class Protein;
class AminoAcid;
class GraphicsProteinHelixItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsProteinHelixItem : public GraphicsItem
{
    public:
        // enumerations
        enum DisplayType{
            Cylinder,
            Ribbon
        };

        // construction and destruction
        GraphicsProteinHelixItem(const QList<const AminoAcid *> &residues);
        ~GraphicsProteinHelixItem();

        // properties
        void setDisplayType(DisplayType type);
        DisplayType displayType() const;

        // drawing
        virtual void paint(GraphicsPainter *painter);

    private:
        GraphicsProteinHelixItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPROTEINHELIXITEM_H
