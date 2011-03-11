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

#ifndef CHEMKIT_GRAPHICSNUCLEICACIDITEM_H
#define CHEMKIT_GRAPHICSNUCLEICACIDITEM_H

#include "graphics.h"

#include "graphicsitem.h"

namespace chemkit {

class Polymer;
class GraphicsNucleicAcidItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsNucleicAcidItem : public GraphicsItem
{
    public:
        // construction and destruction
        GraphicsNucleicAcidItem(const Polymer *polymer = 0);
        ~GraphicsNucleicAcidItem();

        // properties
        void setPolymer(const Polymer *polymer);
        const Polymer* polymer() const;

        // painting
        virtual void paint(GraphicsPainter *painter);

    private:
        GraphicsNucleicAcidItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSNUCLEICACIDITEM_H
