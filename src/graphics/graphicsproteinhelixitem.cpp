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

#include "graphicsproteinhelixitem.h"

#include <chemkit/atom.h>
#include <chemkit/aminoacid.h>

#include "point3g.h"
#include "graphicspainter.h"

namespace chemkit {

// === GraphicsProteinHelixItemPrivate ===================================== //
class GraphicsProteinHelixItemPrivate
{
    public:
        QList<AminoAcid *> residues;
        GraphicsProteinHelixItem::DisplayType displayType;
};

// === GraphicsProteinHelixItem ============================================ //
/// \class GraphicsProteinHelixItem graphicsproteinhelixitem.h chemkit/graphicsproteinhelixitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsProteinHelixItem class visually represents a
///        protein helix.
///
/// \see GraphicsProteinItem

/// \enum GraphicsProteinHelixItem::DisplayType
/// Provides names for the different display types:
///     - \c Cylinder
///     - \c Ribbon

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein helix item object to display \p residues.
GraphicsProteinHelixItem::GraphicsProteinHelixItem(const QList<AminoAcid *> &residues)
    : GraphicsItem(ProteinHelixItem),
      d(new GraphicsProteinHelixItemPrivate)
{
    d->residues = residues;
    d->displayType = Cylinder;
}

/// Destroys the helix item object.
GraphicsProteinHelixItem::~GraphicsProteinHelixItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the display type to \p type.
void GraphicsProteinHelixItem::setDisplayType(DisplayType type)
{
    d->displayType = type;
}

/// Returns the current display type.
GraphicsProteinHelixItem::DisplayType GraphicsProteinHelixItem::displayType() const
{
    return d->displayType;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsProteinHelixItem::paint(GraphicsPainter *painter)
{
    if(d->residues.size() < 2)
        return;

    if(d->displayType == Cylinder){
        if(!d->residues.first()->alphaCarbon() || !d->residues.last()->alphaCarbon())
            return;

        painter->setColor(Qt::red);

        float radius = 1.5;
        Point3g a = d->residues.first()->alphaCarbon()->position();
        Point3g b = d->residues.last()->alphaCarbon()->position();
        painter->drawCylinder(a, b, radius);

        painter->drawCircle(a, radius, (a - b).normalized());
        painter->drawCircle(b, radius, (b - a).normalized());
    }
}

} // end chemkit namespace
