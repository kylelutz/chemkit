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

#include "graphicsproteinitem.h"

#include "graphicsscene.h"
#include "graphicspainter.h"

#include <chemkit/atom.h>
#include <chemkit/foreach.h>
#include <chemkit/polymer.h>
#include <chemkit/residue.h>
#include <chemkit/aminoacid.h>
#include <chemkit/polymerchain.h>

namespace chemkit {

// === GraphicsProteinItemPrivate ========================================== //
class GraphicsProteinItemPrivate
{
public:
    const Polymer *polymer;
    bool secondaryStructureVisible;
    float coilRadius;
    QList<GraphicsProteinCoilItem *> coilItems;
    QList<GraphicsProteinHelixItem *> helixItems;
    QList<GraphicsProteinSheetItem *> sheetItems;
};

// === GraphicsProteinItem ================================================= //
/// \class GraphicsProteinItem graphicsproteinitem.h chemkit/graphicsproteinitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsProteinItem class visually represents a
///        protein polymer.
///
/// GraphicsProteinItem objects manage the following graphics items
/// which display each type of protein secondary structures.
///     - GraphicsProteinCoilItem: Displays random-coil structures
///       using a spline tube.
///     - GraphicsProteinHelixItem: Displays alpha-helicies using
///       cylinders.
///     - GraphicsProteinSheetItem: displays beta-sheets using flat
///       sheets.
///
/// The following image shows a protein item:
/// \image html protein-item.png

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein item to display \p polymer.
GraphicsProteinItem::GraphicsProteinItem(const Polymer *polymer)
    : GraphicsItem(ProteinItem),
      d(new GraphicsProteinItemPrivate)
{
    d->polymer = 0;
    d->secondaryStructureVisible = true;

    setPolymer(polymer);
}

/// Destoys the protein item object.
GraphicsProteinItem::~GraphicsProteinItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the polymer for the protein item to \p polymer.
void GraphicsProteinItem::setPolymer(const Polymer *polymer)
{
    d->polymer = polymer;

    qDeleteAll(d->coilItems);
    d->coilItems.clear();
    qDeleteAll(d->sheetItems);
    d->sheetItems.clear();
    qDeleteAll(d->helixItems);
    d->helixItems.clear();

    if(polymer){
        foreach(const PolymerChain *chain, polymer->chains()){
            if(chain->isEmpty()){
                continue;
            }

            // ensure that the chain contains only amino acids
            bool onlyAminoAcids = true;

            foreach(const Residue *residue, chain->residues()){
                if(residue->residueType() != Residue::AminoAcidResidue){
                    onlyAminoAcids = false;
                    break;
                }
            }

            if(!onlyAminoAcids){
                continue;
            }

            QList<AminoAcid *> residues;
            AminoAcid::Conformation conformation = static_cast<AminoAcid *>(chain->residue(0))->conformation();

            for(size_t i = 0; i < chain->size(); i++){
                AminoAcid *residue = static_cast<AminoAcid *>(chain->residue(i));

                if(residue->conformation() != conformation){
                    if(conformation == AminoAcid::Coil){
                        residues.append(residue);

                        int index = chain->indexOf(residues.first());
                        if(index > 0){
                            residues.prepend(static_cast<AminoAcid *>(chain->residue(index-1)));
                        }

                        GraphicsProteinCoilItem *item = new GraphicsProteinCoilItem(residues);
                        if(scene()){
                            scene()->addItem(item);
                        }

                        d->coilItems.append(item);
                    }
                    else if(conformation == AminoAcid::AlphaHelix){
                        GraphicsProteinHelixItem *item = new GraphicsProteinHelixItem(residues);
                        if(scene()){
                            scene()->addItem(item);
                        }

                        d->helixItems.append(item);
                    }
                    else if(conformation == AminoAcid::BetaSheet){
                        GraphicsProteinSheetItem *item = new GraphicsProteinSheetItem(residues);
                        if(scene()){
                            scene()->addItem(item);
                        }

                        d->sheetItems.append(item);
                    }

                    residues.clear();
                    conformation = residue->conformation();
                }

                residues.append(residue);
            }
        }
    }
}

/// Returns the polymer for the protein item.
const Polymer* GraphicsProteinItem::polymer() const
{
    return d->polymer;
}

/// Sets whether or not the protein's secondary structure is visible.
/// If set to \c false the entire protein will be displayed as if it
/// were all a random-coil structure. The default value is \c true.
void GraphicsProteinItem::setSecondaryStructureVisible(bool visible)
{
    d->secondaryStructureVisible = visible;

    // set the visiblity for all the items
    foreach(GraphicsProteinCoilItem *item, d->coilItems){
        item->setVisible(visible);
    }
    foreach(GraphicsProteinHelixItem *item, d->helixItems){
        item->setVisible(visible);
    }
    foreach(GraphicsProteinSheetItem *item, d->sheetItems){
        item->setVisible(visible);
    }
}

/// Returns \c true if the protein's secondary structure is being
/// displayed.
bool GraphicsProteinItem::secondaryStructureVisible() const
{
    return d->secondaryStructureVisible;
}

// --- Painting ------------------------------------------------------------ //
void GraphicsProteinItem::paint(GraphicsPainter *painter)
{
    if(!d->polymer || d->polymer->chainCount() < 1){
        return;
    }

    // draw a single spline through all of the residues' alpha-carbons for each chain
    if(!d->secondaryStructureVisible){
        float radius = 0.35f;

        foreach(const PolymerChain *chain, d->polymer->chains()){
            QList<Point3f> trace;

            foreach(const Residue *residue, chain->residues()){
                const AminoAcid *aminoAcid = static_cast<const AminoAcid *>(residue);

                if(aminoAcid->alphaCarbon() != 0){
                    trace.append(aminoAcid->alphaCarbon()->position().cast<float>());
                }
            }

            painter->setColor(Qt::green);
            painter->drawSpline(trace, radius, 3);
        }
    }
}

} // end chemkit namespace
