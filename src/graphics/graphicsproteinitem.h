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
    void setCoilRadius(float radius);
    float coilRadius() const;
    void setHelixDisplayType(GraphicsProteinHelixItem::DisplayType type);
    GraphicsProteinHelixItem::DisplayType helixDisplayType() const;

    // painting
    virtual void paint(GraphicsPainter *painter);

private:
    GraphicsProteinItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPROTEINITEM_H
