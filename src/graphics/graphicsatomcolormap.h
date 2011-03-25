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

#ifndef CHEMKIT_GRAPHICSATOMCOLORMAP_H
#define CHEMKIT_GRAPHICSATOMCOLORMAP_H

#include "graphics.h"

#include <QColor>

namespace chemkit {

class Atom;
class Element;
class GraphicsAtomColorMapPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsAtomColorMap
{
    public:
        // enumerations
        enum ColorScheme {
            DefaultColorScheme,
            RasmolColorScheme,
            PymolColorScheme,
            JmolColorScheme
        };

        // construction and destruction
        GraphicsAtomColorMap();
        GraphicsAtomColorMap(ColorScheme scheme);
        GraphicsAtomColorMap(const GraphicsAtomColorMap &colorMap);
        virtual ~GraphicsAtomColorMap();

        // colors
        void setColor(const Element &element, const QColor &color);
        virtual QColor color(const Element &element) const;
        virtual QColor color(const Atom *atom) const;
        void setDefaultColor(const QColor &color);
        QColor defaultColor() const;
        void setColorScheme(ColorScheme scheme);

    private:
        GraphicsAtomColorMapPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSATOMCOLORMAP_H
