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

#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/polymer.h>
#include <chemkit/nucleotide.h>
#include <chemkit/polymerchain.h>

namespace chemkit {

// === GraphicsNucleicAcidItemPrivate ====================================== //
class GraphicsNucleicAcidItemPrivate
{
    public:
        const Polymer *polymer;
};

// === GraphicsNucleicAcidItem ============================================= //
/// \class GraphicsNucleicAcidItem graphicsnucleicaciditem.h chemkit/graphicsnucleicaciditem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsNucleicAcidItem represents a nucleic acid.
///
/// The GraphicsNucleicAcid item displays a nucleic acid Polymer.
///
/// \image html nucleic-acid-item.png

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new nucleic acid item to display \p polymer.
GraphicsNucleicAcidItem::GraphicsNucleicAcidItem(const Polymer *polymer)
    : GraphicsItem(NucleicAcidItem),
      d(new GraphicsNucleicAcidItemPrivate)
{
    d->polymer = polymer;
}

/// Destroys the nucleic acid item.
GraphicsNucleicAcidItem::~GraphicsNucleicAcidItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the polymer for the item to display to \p polymer.
void GraphicsNucleicAcidItem::setPolymer(const Polymer *polymer)
{
    d->polymer = polymer;
}

/// Returns the polymer for the item.
const Polymer* GraphicsNucleicAcidItem::polymer() const
{
    return d->polymer;
}

// --- Painting ------------------------------------------------------------ //
void GraphicsNucleicAcidItem::paint(GraphicsPainter *painter)
{
    if(!d->polymer){
        return;
    }

    foreach(const PolymerChain *chain, d->polymer->chains()){
        // ensure the chain contains only nucleotides
        bool isOnlyNucleotides = true;

        foreach(const Residue *residue, chain->residues()){
            if(residue->residueType() != Residue::NucleotideResidue){
                isOnlyNucleotides = false;
                break;
            }
        }

        if(!isOnlyNucleotides){
            continue;
        }

        // build a list of points for the spline
        QList<Point3f> trace;
        float radius = 0.6f;
        painter->setColor(Qt::cyan);

        foreach(const Residue *residue, chain->residues()){
            // add phosphorus position to trace
            const Atom *phosphorus = residue->atom("P");
            if(phosphorus){
                trace.append(phosphorus->position());

                // draw ladder
                float ladderRadius = 0.4f;
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
