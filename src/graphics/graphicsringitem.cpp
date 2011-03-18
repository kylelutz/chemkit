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

#include "graphicsringitem.h"

#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/ring.h>

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

    Point3g center = ringCenter(d->ring);
    QList<Atom *> atoms = d->ring->atoms();

    for(int i = 0; i < atoms.size(); i++){
        const Atom *a = atoms[i];
        const Atom *b = atoms[(i + 1) % d->ring->size()];

        painter->drawTriangle(a->position(), b->position(), center);
        painter->drawTriangle(b->position(), a->position(), center);
    }
}

// --- Internal Methods ---------------------------------------------------- //
Point3g GraphicsRingItem::ringCenter(const Ring *ring) const
{
    GraphicsFloat sx = 0;
    GraphicsFloat sy = 0;
    GraphicsFloat sz = 0;

    foreach(const Atom *atom, ring->atoms()){
        sx += atom->x();
        sy += atom->y();
        sz += atom->z();
    }

    int n = ring->atomCount();

    return Point3g(sx/n, sy/n, sz/n);
}

} // end chemkit namespace
