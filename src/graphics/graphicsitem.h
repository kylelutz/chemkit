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

#ifndef CHEMKIT_GRAPHICSITEM_H
#define CHEMKIT_GRAPHICSITEM_H

#include "graphics.h"

#include "graphicstransform.h"

namespace chemkit {

class GraphicsRay;
class GraphicsScene;
class GraphicsPainter;
class GraphicsMaterial;
class GraphicsItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsItem
{
    public:
        // enumerations
        enum ItemType {
            GenericItem,
            AtomItem,
            BondItem,
            MoleculeItem,
            ProteinItem,
            ProteinCoilItem,
            ProteinHelixItem,
            ProteinSheetItem,
            NucleicAcidItem,
            CustomItem = 256
        };

        enum ItemChange {
            ItemOpacityChanged,
            ItemSceneChanged,
            ItemVisiblityChanged
        };

        // construction and destruction
        GraphicsItem(int type = GenericItem);
        virtual ~GraphicsItem();

        // properties
        int type() const;
        void setVisible(bool visible);
        bool isVisible() const;
        void show();
        void hide();
        GraphicsScene* scene() const;
        void setOpacity(GraphicsFloat opacity);
        GraphicsFloat opacity() const;
        bool isOpaque() const;
        bool isTransparent() const;
        bool isTranslucent() const;

        // material
        void setMaterial(GraphicsMaterial *material);
        GraphicsMaterial* material() const;

        // geometry
        void setTransform(const GraphicsTransform &transform);
        GraphicsTransform transform() const;
        void translate(const Vector3g &vector);
        void translate(GraphicsFloat x, GraphicsFloat y, GraphicsFloat z);
        void rotate(const Vector3g &axis, const GraphicsFloat angle);
        virtual bool intersects(const GraphicsRay &ray, GraphicsFloat *distance = 0) const;

        // drawing
        virtual void paint(GraphicsPainter *painter);
        void update();

    protected:
        // events
        virtual void itemChanged(ItemChange change);

    private:
        // internal methods
        void setScene(GraphicsScene *scene);

        friend class GraphicsScene;

    private:
        GraphicsItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSITEM_H
