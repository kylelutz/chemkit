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

#include "graphicsproteincoilitem.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/foreach.h>
#include <chemkit/aminoacid.h>

#include "graphicspainter.h"

namespace chemkit {

// === GraphicsProteinCoilItemPrivate ====================================== //
class GraphicsProteinCoilItemPrivate
{
public:
    int curveDegree;
    float radius;
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
void GraphicsProteinCoilItem::setRadius(float radius)
{
    d->radius = radius;
}

/// Returns the radius of the coil.
float GraphicsProteinCoilItem::radius() const
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
    QList<Point3f> trace;
    foreach(const AminoAcid *residue, d->residues){
        const Atom *alphaCarbon = residue->alphaCarbon();
        if(alphaCarbon){
            trace.append(alphaCarbon->position().cast<float>());
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
