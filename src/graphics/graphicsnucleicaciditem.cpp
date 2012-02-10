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

#include "graphicsnucleicaciditem.h"

#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/foreach.h>
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
                trace.append(phosphorus->position().cast<float>());

                // draw ladder
                float ladderRadius = 0.4f;
                const Atom *centerAtom = residue->atom("C2");
                if(centerAtom){
                    painter->drawCylinder(phosphorus->position().cast<float>(),
                                          centerAtom->position().cast<float>(),
                                          ladderRadius);

                    painter->drawSphere(centerAtom->position().cast<float>(),
                                        ladderRadius);
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
