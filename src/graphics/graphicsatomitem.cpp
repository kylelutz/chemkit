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

#include "graphicsatomitem.h"

#include "graphicsray.h"
#include "graphicssphere.h"
#include "graphicspainter.h"
#include "graphicsmoleculeitem.h"
#include "graphicsvertexbuffer.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>

namespace chemkit {

class GraphicsAtomItemPrivate
{
    public:
        const Atom *atom;
        GraphicsSphere sphere;
        QColor color;
        GraphicsVertexBuffer *vertexBuffer;
};

// === GraphicsAtomItem ==================================================== //
/// \class GraphicsAtomItem graphicsatomitem.h chemkit/graphicsatomitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsAtomItem class displays a single atom.
///
/// The GraphicsAtomItem class displays a single atom using a sphere.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new atom item that displays \p atom with a sphere of
/// the given \p radius.
GraphicsAtomItem::GraphicsAtomItem(const Atom *atom, float radius)
    : GraphicsItem(AtomItem),
      d(new GraphicsAtomItemPrivate)
{
    d->vertexBuffer = 0;
    d->atom = atom;
    d->sphere = GraphicsSphere(radius);

    if(atom){
        translate(atom->position().cast<float>());
    }
}

/// Destroys the atom item.
GraphicsAtomItem::~GraphicsAtomItem()
{
    delete d->vertexBuffer;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the atom the item displays to \p atom.
void GraphicsAtomItem::setAtom(const Atom *atom)
{
    d->atom = atom;

    if(atom){
        setTransform(GraphicsTransform::translation(atom->position().cast<float>()));
    }
}

/// Returns the atom the item displays.
const Atom* GraphicsAtomItem::atom() const
{
    return d->atom;
}

/// Sets the radius of the sphere used to display the atom.
void GraphicsAtomItem::setRadius(float radius)
{
    d->sphere.setRadius(radius);

    delete d->vertexBuffer;
    d->vertexBuffer = 0;
}

/// Returns the radius of the sphere used to display the atom.
float GraphicsAtomItem::radius() const
{
    return d->sphere.radius();
}

/// Sets the color for the atom item to \p color.
void GraphicsAtomItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color for the atom item.
QColor GraphicsAtomItem::color() const
{
    return d->color;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns \c true if the item intersects \p ray.
bool GraphicsAtomItem::intersects(const GraphicsRay &ray, float *distance) const
{
    return ray.intersectsSphere(d->atom->position().cast<float>(), radius(), distance);
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsAtomItem::paint(GraphicsPainter *painter)
{
    if(!d->vertexBuffer){
        d->vertexBuffer = d->sphere.tesselate();
    }

    painter->setColor(d->color);
    painter->draw(d->vertexBuffer);
}

} // end chemkit namespace
