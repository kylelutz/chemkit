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

#include "graphicsnucleicaciditem.h"

#include "graphicspoint.h"
#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/nucleotide.h>
#include <chemkit/nucleicacid.h>
#include <chemkit/nucleicacidchain.h>

namespace chemkit {

// === GraphicsNucleicAcidItemPrivate ====================================== //
class GraphicsNucleicAcidItemPrivate
{
    public:
        const NucleicAcid *nucleicAcid;
};

// === GraphicsNucleicAcidItem ============================================= //
/// \class GraphicsNucleicAcidItem graphicsnucleicaciditem.h chemkit/graphicsnucleicaciditem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsNucleicAcidItem represents a nucleic acid.
///
/// The GraphicsNucleicAcid item displays a NucleicAcid object.
///
/// \image html nucleic-acid-item.png

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new nucleic acid item to display \p nucleicAcid.
GraphicsNucleicAcidItem::GraphicsNucleicAcidItem(const NucleicAcid *nucleicAcid)
    : GraphicsItem(NucleicAcidItem),
      d(new GraphicsNucleicAcidItemPrivate)
{
    d->nucleicAcid = nucleicAcid;
}

/// Destroys the nucleic acid item.
GraphicsNucleicAcidItem::~GraphicsNucleicAcidItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the nucleic acid to display.
void GraphicsNucleicAcidItem::setNucleicAcid(const NucleicAcid *nucleicAcid)
{
    d->nucleicAcid = nucleicAcid;
}

/// Returns the current nucleic acid item or \c 0 if no nucleic aicd
/// is being displayed.
const NucleicAcid* GraphicsNucleicAcidItem::nucleicAcid() const
{
    return d->nucleicAcid;
}

// --- Painting ------------------------------------------------------------ //
void GraphicsNucleicAcidItem::paint(GraphicsPainter *painter)
{
    if(!d->nucleicAcid){
        return;
    }

    foreach(const NucleicAcidChain *chain, d->nucleicAcid->chains()){
        GraphicsFloat radius = 0.6f;

        QList<GraphicsPoint> trace;
        painter->setColor(Qt::cyan);

        foreach(const Nucleotide *residue, chain->residues()){
            // add phosphorus position to trace
            const Atom *phosphorus = residue->atom("P");
            if(phosphorus){
                trace.append(phosphorus->position());

                // draw ladder
                GraphicsFloat ladderRadius = 0.4f;
                const Atom *centerAtom = residue->atom("C2");
                if(centerAtom){
                    painter->drawCylinder(phosphorus->position(), centerAtom->position(), ladderRadius);
                    painter->drawSphere(centerAtom->position(), ladderRadius);
                }
            }
        }

        if(trace.size() > 2){
            painter->setColor(Qt::blue);
            painter->drawSpline(trace, radius, 3);

            painter->drawSphere(trace.first(), radius);
            painter->drawSphere(trace.last(), radius);
        }
    }
}

} // end chemkit namespace
