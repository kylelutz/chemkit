/******************************************************************************
**
** Copyright (C) 2011 Wang Wenlin <sopl.wang@gmail.com>
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

#ifndef CHEMKIT_GRAPHICSPYMOLSURFACEITEM_H
#define CHEMKIT_GRAPHICSPYMOLSURFACEITEM_H

#include "graphics.h"

#include <boost/shared_ptr.hpp>

#include "graphicsitem.h"

#include <chemkit/atomcolormap.h>

namespace chemkit {

class Molecule;
class GraphicsPymolSurfaceItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsPymolSurfaceItem : public GraphicsItem
{
public:
    // enumerations
    enum SurfaceQuality {
        SurfaceQualityMoreMiserable = -4,
        SurfaceQualityMiserable = -3,
        SurfaceQualityMorePoor = -2,
        SurfaceQualityPoor = -1,
        SurfaceQualityNormal = 0,
        SurfaceQualityGood = 1,
        SurfaceQualityNearPerfect = 2,
        SurfaceQualityPerfect = 3,
        SurfaceQualityImpractical = 4
    };

    enum SurfaceType {
        SurfaceTypeSolid = 0,
        SurfaceTypeDots = 1,
        SurfaceTypeTriangles = 2,
        SurfaceTypeType3 = 3,
        SurfaceTypeType4 = 4,
        SurfaceTypeType5 = 5,
        SurfaceTypeType6 = 6
    };

    enum SolventType {
        SolventTypeExcluded = 0,
        SolventTypeAccessible = 1
    };

    enum ColorMode {
        SolidColor,
        AtomColor
    };

    // construction and destruction
    GraphicsPymolSurfaceItem(const Molecule *molecule = 0, SolventType solventType = SolventTypeExcluded);
    ~GraphicsPymolSurfaceItem();

    // properties
    void setMolecule(const Molecule *molecule);
    const Molecule* molecule() const;
    void setQuality(SurfaceQuality quality);
    SurfaceQuality quality() const;
    void setSurfaceType(SurfaceType type);
    SurfaceType surfaceType() const;
    void setSolventType(SolventType solventType);
    SolventType solventType() const;
    void setProbeRadius(Real radius);
    Real probeRadius() const;
    void setColorMode(ColorMode mode);
    ColorMode colorMode() const;
    void setColor(const QColor &color);
    QColor color() const;
    void setColorMap(const boost::shared_ptr<AtomColorMap> &colorMap);
    boost::shared_ptr<AtomColorMap> colorMap() const;

    // drawing
    virtual void paint(GraphicsPainter *painter);

private:
    // internal methods
    void itemChanged(ItemChange change);
    void setCalculated(bool calculated);
    Real maxVdwRadius();

private:
    GraphicsPymolSurfaceItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPYMOLSURFACEITEM_H
