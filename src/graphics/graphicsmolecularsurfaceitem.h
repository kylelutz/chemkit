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

#ifndef CHEMKIT_GRAPHICSMOLECULARSURFACEITEM_H
#define CHEMKIT_GRAPHICSMOLECULARSURFACEITEM_H

#include "graphics.h"

#include "graphicsitem.h"
#include "graphicsatomcolormap.h"

#include <chemkit/molecularsurface.h>

namespace chemkit {

class GraphicsMolecularSurfaceItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsMolecularSurfaceItem : public GraphicsItem
{
    public:
        // enumerations
        enum ColorMode {
            SolidColor,
            AtomColor
        };

        // construction and destruction
        GraphicsMolecularSurfaceItem(const Molecule *molecule = 0);
        GraphicsMolecularSurfaceItem(const MolecularSurface *surface);
        ~GraphicsMolecularSurfaceItem();

        // properties
        void setSurface(const MolecularSurface *surface);
        const MolecularSurface* surface() const;
        void setMolecule(const Molecule *molecule);
        const Molecule* molecule() const;
        void setSurfaceType(MolecularSurface::SurfaceType type);
        MolecularSurface::SurfaceType surfaceType() const;
        void setProbeRadius(float radius);
        float probeRadius() const;
        void setColor(const QColor &color);
        QColor color() const;
        void setColorMode(ColorMode mode);
        ColorMode colorMode() const;
        void setAtomColorMap(GraphicsAtomColorMap *colorMap);
        GraphicsAtomColorMap* colorMap() const;

    private:
        void itemChanged(ItemChange change);
        void recalculate();

    private:
        GraphicsMolecularSurfaceItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSMOLECULARSURFACEITEM_H
