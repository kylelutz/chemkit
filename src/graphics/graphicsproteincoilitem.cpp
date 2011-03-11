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

#include "graphicsproteincoilitem.h"

#include <chemkit/atom.h>
#include <chemkit/aminoacid.h>

#include "graphicspoint.h"
#include "graphicspainter.h"

namespace chemkit {

// === GraphicsProteinCoilItemPrivate ====================================== //
class GraphicsProteinCoilItemPrivate
{
    public:
        int curveDegree;
        GraphicsFloat radius;
        QList<AminoAcid *> residues;
};

// === GraphicsProteinCoilItem ============================================= //
/// \class GraphicsProteinCoilItem graphicsproteincoilitem.h chemkit/graphicsproteincoilitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsProteinCoilItem class visually represents a
///        protein coil.
///
/// \see GraphicsProteinItem

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein coil item to display \p residues.
GraphicsProteinCoilItem::GraphicsProteinCoilItem(const QList<AminoAcid *> &residues)
    : GraphicsItem(ProteinCoilItem),
      d(new GraphicsProteinCoilItemPrivate)
{
    d->residues = residues;
    d->radius = 0.35f;
    d->curveDegree = 2;
}

/// Destroys the protein coil item object.
GraphicsProteinCoilItem::~GraphicsProteinCoilItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the radius of the coil to \p radius.
void GraphicsProteinCoilItem::setRadius(GraphicsFloat radius)
{
    d->radius = radius;
}

/// Returns the radius of the coil.
GraphicsFloat GraphicsProteinCoilItem::radius() const
{
    return d->radius;
}

/// Sets the curve degree to \p degree.
void GraphicsProteinCoilItem::setCurveDegree(int degree)
{
    d->curveDegree = degree;
}

/// Returns the curve degree.
int GraphicsProteinCoilItem::curveDegree() const
{
    return d->curveDegree;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsProteinCoilItem::paint(GraphicsPainter *painter)
{
    // build list of positions of all the alpha carbons
    QList<GraphicsPoint> trace;
    foreach(const AminoAcid *residue, d->residues){
        const Atom *alphaCarbon = residue->alphaCarbon();
        if(alphaCarbon){
            trace.append(alphaCarbon->position());
        }
    }

    // nothing to do if the trace contains less that 2 points
    if(trace.size() < 2){
        return;
    }

    // draw the spline
    painter->setColor(Qt::green);
    painter->drawSpline(trace, d->radius, d->curveDegree + 1);
}

} // end chemkit namespace
