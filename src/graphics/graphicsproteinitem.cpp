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

#include "graphicsproteinitem.h"

#include "graphicsscene.h"
#include "graphicspainter.h"

#include <chemkit/protein.h>

namespace chemkit {

// === GraphicsProteinItemPrivate ========================================== //
class GraphicsProteinItemPrivate
{
    public:
        const Protein *protein;
        bool secondaryStructureVisible;
        GraphicsFloat coilRadius;
        QList<GraphicsProteinCoilItem *> coilItems;
        QList<GraphicsProteinHelixItem *> helixItems;
        QList<GraphicsProteinSheetItem *> sheetItems;
};

// === GraphicsProteinItem ================================================= //
/// \class GraphicsProteinItem graphicsproteinitem.h chemkit/graphicsproteinitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsProteinItem class visually represents a
///        protein.
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
/// Creates a new protein item to display \p protein.
GraphicsProteinItem::GraphicsProteinItem(const Protein *protein)
    : GraphicsItem(ProteinItem),
      d(new GraphicsProteinItemPrivate)
{
    d->protein = protein;
    d->secondaryStructureVisible = true;
}

/// Destoys the protein item object.
GraphicsProteinItem::~GraphicsProteinItem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the protein to display.
void GraphicsProteinItem::setProtein(const Protein *protein)
{
    d->protein = protein;

    qDeleteAll(d->coilItems);
    d->coilItems.clear();
    qDeleteAll(d->sheetItems);
    d->sheetItems.clear();
    qDeleteAll(d->helixItems);
    d->helixItems.clear();

    if(protein){
        foreach(const ProteinChain *chain, protein->chains()){
            if(chain->residueCount() < 1){
                continue;
            }

            QList<AminoAcid *> residues;
            AminoAcid::Conformation conformation = chain->residue(0)->conformation();

            for(int i = 0; i < chain->size(); i++){
                AminoAcid *residue = chain->residue(i);

                if(residue->conformation() != conformation){
                    if(conformation == AminoAcid::Coil){
                        residues.append(residue);

                        int index = chain->residues().indexOf(residues.first());
                        if(index > 0)
                            residues.prepend(chain->residue(index-1));

                        GraphicsProteinCoilItem *item = new GraphicsProteinCoilItem(residues);
                        if(scene())
                            scene()->addItem(item);
                        d->coilItems.append(item);
                    }
                    else if(conformation == AminoAcid::AlphaHelix){
                        GraphicsProteinHelixItem *item = new GraphicsProteinHelixItem(residues);
                        if(scene())
                            scene()->addItem(item);
                        d->helixItems.append(item);
                    }
                    else if(conformation == AminoAcid::BetaSheet){
                        GraphicsProteinSheetItem *item = new GraphicsProteinSheetItem(residues);
                        if(scene())
                            scene()->addItem(item);
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

/// Returns the protein that is being displayed.
const Protein* GraphicsProteinItem::protein() const
{
    return d->protein;
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
    if(!d->protein || d->protein->chainCount() < 1){
        return;
    }

    // draw a single spline through all of the residues' alpha-carbons for each chain
    if(!d->secondaryStructureVisible){
        GraphicsFloat radius = 0.35f;

        foreach(const ProteinChain *chain, d->protein->chains()){
            QList<GraphicsPoint> trace;

            foreach(const AminoAcid *residue, chain->residues()){
                if(residue->alphaCarbon() != 0){
                    trace.append(residue->alphaCarbon()->position());
                }
            }

            painter->setColor(Qt::green);
            painter->drawSpline(trace, radius, 3);
        }
    }
}

} // end chemkit namespace
