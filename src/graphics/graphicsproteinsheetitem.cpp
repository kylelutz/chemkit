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

#include "graphicsproteinsheetitem.h"

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/aminoacid.h>

#include "graphicspainter.h"

namespace chemkit {

// === GraphicsProteinSheetItemPrivate ===================================== //
class GraphicsProteinSheetItemPrivate
{
public:
    QList<AminoAcid *> residues;
};

// === GraphicsProteinSheetItem ============================================ //
/// \class GraphicsProteinSheetItem graphicsproteinsheetitem.h chemkit/graphicsproteinsheetitem.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsProteinSheetItem class visually represents a
///        protein beta-sheet.
///
/// \see GraphicsProteinItem

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein sheet item to display \p residues.
GraphicsProteinSheetItem::GraphicsProteinSheetItem(const QList<AminoAcid *> &residues)
    : GraphicsItem(ProteinSheetItem),
      d(new GraphicsProteinSheetItemPrivate)
{
    d->residues = residues;
}

/// Destroys the protein sheet item.
GraphicsProteinSheetItem::~GraphicsProteinSheetItem()
{
    delete d;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsProteinSheetItem::paint(GraphicsPainter *painter)
{
    if(d->residues.size() < 2){
        return;
    }

    QVector<Point3f> trace(d->residues.size());
    QVector<Vector3f> normals(d->residues.size());

    for(int i = 0; i < d->residues.size(); i++){
        const AminoAcid *residue = d->residues[i];
        trace[i] = residue->alphaCarbon()->position().cast<float>();

        Vector3f normal = residue->peptidePlaneNormal().cast<float>();

        // flip every other normal
        if(i & 1){
            normal = -normal;
        }

        normals[i] = normal.normalized();
    }

    // surface settings
    float width = 2.0;
    float height = width / 4.0;

    int degree = 3;
    int order = degree + 1;

    QVector<Point3f> controlPoints(trace.size() * 5);

    for(int i = 0; i < trace.size(); i++){
        const Point3f &point = trace[i];
        const Vector3f &normal = normals[i];

        Vector3f axis;
        if(i == 0)
            axis = trace[i+1] - trace[i];
        else
            axis = trace[i] - trace[i-1];

        Vector3f right = axis.cross(normal).normalized();

        float halfHeight = height / 2.0;
        float halfWidth = width / 2.0;

        // four points for the rectangle
        // counter clockwise winding (topLeft -> bottomLeft -> bottomRight -> topRight)

        // top left
        Point3f topLeft = point;
        topLeft += normal * halfHeight;
        topLeft += right * -halfWidth;
        controlPoints[i*5+0] = topLeft;

        // bottom left
        Point3f bottomLeft = point;
        bottomLeft += normal * -halfHeight;
        bottomLeft += right * -halfWidth;
        controlPoints[i*5+1] = bottomLeft;

        // bottom right
        Point3f bottomRight = point;
        bottomRight += normal * -halfHeight;
        bottomRight += right * halfWidth;
        controlPoints[i*5+2] = bottomRight;

        // top right
        Point3f topRight = point;
        topRight += normal * halfHeight;
        topRight += right * halfWidth;
        controlPoints[i*5+3] = topRight;

        // top left (again)
        controlPoints[i*5+4] = topLeft;
    }

    QVector<float> uKnots(trace.size() + degree + 1);
    for(int i = 0; i < uKnots.size(); i++){
        if(i < order)
            uKnots[i] = 0;
        else if(i > (uKnots.size() - 1 - order))
            uKnots[i] = uKnots.size()-2*order+1;
        else
            uKnots[i] = i-degree;
    }

    QVector<float> vKnots(7);
    vKnots[0] = 0;
    vKnots[1] = 0;
    vKnots[2] = 1;
    vKnots[3] = 2;
    vKnots[4] = 3;
    vKnots[5] = 4;
    vKnots[6] = 4;

    // draw the nurbs surface
    painter->setColor(Qt::yellow);
    painter->drawNurbsSurface(controlPoints, uKnots, vKnots, order, 2);

    // draw rectangle at start
    painter->drawRectangle(controlPoints[0],
                           controlPoints[1],
                           controlPoints[2],
                           controlPoints[3]);

    // draw rectangle at end
    painter->drawRectangle(controlPoints[controlPoints.size()-1],
                           controlPoints[controlPoints.size()-2],
                           controlPoints[controlPoints.size()-3],
                           controlPoints[controlPoints.size()-4]);
}

} // end chemkit namespace
