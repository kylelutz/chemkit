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

#include "graphicsringitem.h"

#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/ring.h>
#include <chemkit/foreach.h>

namespace chemkit {

// === GraphicsRingItemPrivate ============================================= //
class GraphicsRingItemPrivate
{
public:
    const Ring *ring;
    QColor color;
};

// === GraphicsRingItem ==================================================== //
/// \class GraphicsRingItem graphicsringitem.h chemkit/graphicsringitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsRingItem class displays a ring.
///
/// The image below shows two ring items in a uridine molecule:
/// \image html ring-item.png

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new ring item to display \p ring.
GraphicsRingItem::GraphicsRingItem(const Ring *ring)
    : GraphicsItem(),
      d(new GraphicsRingItemPrivate)
{
    d->ring = ring;
    d->color = Qt::blue;
}

/// Destroys the ring item object.
GraphicsRingItem::~GraphicsRingItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the ring for the item to display to \p ring.
void GraphicsRingItem::setRing(const Ring *ring)
{
    d->ring = ring;
}

/// Returns the ring that the item is displaying.
const Ring* GraphicsRingItem::ring() const
{
    return d->ring;
}

/// Sets the color of the ring item to \p color.
void GraphicsRingItem::setColor(const QColor &color)
{
    d->color = color;
}

/// Returns the color of the ring item.
QColor GraphicsRingItem::color() const
{
    return d->color;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsRingItem::paint(GraphicsPainter *painter)
{
    if(!d->ring){
        return;
    }

    QColor color = d->color;
    color.setAlphaF(opacity());
    painter->setColor(color);

    Point3f center = ringCenter(d->ring);
    std::vector<Atom *> atoms(d->ring->atoms().begin(),
                              d->ring->atoms().end());

    for(unsigned int i = 0; i < atoms.size(); i++){
        const Atom *a = atoms[i];
        const Atom *b = atoms[(i + 1) % d->ring->size()];

        painter->drawTriangle(a->position().cast<float>(), b->position().cast<float>(), center);
        painter->drawTriangle(b->position().cast<float>(), a->position().cast<float>(), center);
    }
}

// --- Internal Methods ---------------------------------------------------- //
Point3f GraphicsRingItem::ringCenter(const Ring *ring) const
{
    float sx = 0;
    float sy = 0;
    float sz = 0;

    foreach(const Atom *atom, ring->atoms()){
        sx += atom->x();
        sy += atom->y();
        sz += atom->z();
    }

    int n = ring->atomCount();

    return Point3f(sx/n, sy/n, sz/n);
}

} // end chemkit namespace
