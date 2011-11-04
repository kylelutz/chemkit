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

#include "graphicsproteinhelixitem.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/aminoacid.h>

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
        Point3f a = d->residues.first()->alphaCarbon()->position().cast<float>();
        Point3f b = d->residues.last()->alphaCarbon()->position().cast<float>();
        painter->drawCylinder(a, b, radius);

        painter->drawCircle(a, radius, (a - b).normalized());
        painter->drawCircle(b, radius, (b - a).normalized());
    }
}

} // end chemkit namespace
