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

#ifndef CHEMKIT_GRAPHICSMOLECULEITEM_H
#define CHEMKIT_GRAPHICSMOLECULEITEM_H

#include "graphics.h"

#include <boost/shared_ptr.hpp>

#include <chemkit/molecule.h>
#include <chemkit/atomcolormap.h>
#include <chemkit/moleculewatcher.h>

#include "graphicsitem.h"
#include "graphicsatomitem.h"
#include "graphicsbonditem.h"

namespace chemkit {

class GraphicsMoleculeItemPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsMoleculeItem : public GraphicsItem
{
public:
    // enumerations
    enum DisplayType {
        BallAndStick,
        Stick,
        SpaceFilling
    };

    // construction and destruction
    GraphicsMoleculeItem(const Molecule *molecule = 0);
    ~GraphicsMoleculeItem();

    // properties
    void setMolecule(const Molecule *molecule);
    const Molecule* molecule() const;
    void setDisplayType(DisplayType type);
    DisplayType displayType() const;
    void setAtomRadius(float radius);
    float atomRadius() const;
    void setBondRadius(float radius);
    float bondRadius() const;
    void setAtomVisible(const Atom *atom, bool visible);
    bool atomVisible(const Atom *atom) const;
    void setHydrogensVisible(bool visible);
    bool hydrogensVisible() const;
    void setHydrogenScale(float scale);
    float hydrogenScale() const;
    void setBondOrderVisible(bool visible);
    bool bondOrderVisible() const;
    void setAtomColoredBonds(bool atomColoredBonds);
    bool atomColoredBonds() const;
    void setColorMap(const boost::shared_ptr<AtomColorMap> &colorMap);
    boost::shared_ptr<AtomColorMap> colorMap() const;

    // items
    GraphicsAtomItem* atomItem(const Atom *atom) const;
    GraphicsBondItem* bondItem(const Bond *bond) const;

    // painting
    virtual void paint(GraphicsPainter *painter);

protected:
    // events
    void itemChanged(ItemChange change);

private:
    // slots
    void atomAdded(const Atom *atom);
    void atomRemoved(const Atom *atom);
    void atomElementChanged(const Atom *atom);
    void atomPositionChanged(const Atom *atom);
    void bondAdded(const Bond *bond);
    void bondRemoved(const Bond *bond);
    void bondOrderChanged(const Bond *bond);

private:
    GraphicsMoleculeItemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSMOLECULEITEM_H
