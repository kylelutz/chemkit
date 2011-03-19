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

#include "graphicsatomitem.h"

#include "point3g.h"
#include "graphicsray.h"
#include "graphicssphere.h"
#include "graphicspainter.h"
#include "graphicsmoleculeitem.h"
#include "graphicsvertexbuffer.h"

#include <chemkit/atom.h>

namespace chemkit {

class GraphicsAtomItemPrivate
{
    public:
        const Atom *atom;
        GraphicsSphere sphere;
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
        translate(Point3g(atom->position()));
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
        setTransform(GraphicsTransform::translation(Point3g(atom->position())));
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

// --- Geometry ------------------------------------------------------------ //
/// Returns \c true if the item intersects \p ray.
bool GraphicsAtomItem::intersects(const GraphicsRay &ray, float *distance) const
{
    return ray.intersectsSphere(d->atom->position(), radius(), distance);
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsAtomItem::paint(GraphicsPainter *painter)
{
    if(!d->vertexBuffer){
        d->vertexBuffer = d->sphere.tesselate();
    }

    painter->setColor(GraphicsMoleculeItem::atomColor(d->atom));
    painter->draw(d->vertexBuffer);
}

} // end chemkit namespace
