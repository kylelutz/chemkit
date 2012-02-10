/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#ifndef CHEMKIT_GRAPHICSITEM_H
#define CHEMKIT_GRAPHICSITEM_H

#include "graphics.h"

#include "graphicstransform.h"
#include "graphicsboundingbox.h"

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
    void setOpacity(float opacity);
    float opacity() const;
    bool isOpaque() const;
    bool isTransparent() const;
    bool isTranslucent() const;

    // material
    void setMaterial(GraphicsMaterial *material);
    GraphicsMaterial* material() const;

    // geometry
    void setTransform(const GraphicsTransform &transform);
    GraphicsTransform transform() const;
    void translate(const Vector3f &vector);
    void translate(float x, float y, float z);
    void rotate(const Vector3f &axis, const float angle);
    virtual GraphicsBoundingBox boundingBox() const;
    virtual bool intersects(const GraphicsRay &ray, float *distance = 0) const;

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
