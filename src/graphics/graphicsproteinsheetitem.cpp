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

#include "graphicsproteinsheetitem.h"

#include <chemkit/atom.h>
#include <chemkit/aminoacid.h>

#include "point3g.h"
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

    QVector<Point3g> trace(d->residues.size());
    QVector<GraphicsVector> normals(d->residues.size());

    for(int i = 0; i < d->residues.size(); i++){
        const AminoAcid *residue = d->residues[i];
        trace[i] = residue->alphaCarbon()->position();

        GraphicsVector normal = GraphicsVector(residue->peptidePlaneNormal());

        // flip every other normal
        if(i & 1){
            normal = -normal;
        }

        normals[i] = normal.normalized();
    }

    // surface settings
    GraphicsFloat width = 2.0;
    GraphicsFloat height = width / 4.0;

    int degree = 3;
    int order = degree + 1;

    QVector<Point3g> controlPoints(trace.size() * 5);

    for(int i = 0; i < trace.size(); i++){
        const Point3g &point = trace[i];
        const GraphicsVector &normal = normals[i];

        GraphicsVector axis;
        if(i == 0)
            axis = trace[i+1] - trace[i];
        else
            axis = trace[i] - trace[i-1];

        GraphicsVector right = axis.cross(normal);

        GraphicsFloat halfHeight = height / 2.0;
        GraphicsFloat halfWidth = width / 2.0;

        // four points for the rectangle
        // counter clockwise winding (topLeft -> bottomLeft -> bottomRight -> topRight)

        // top left
        Point3g topLeft = point;
        topLeft.moveBy(halfHeight, normal);
        topLeft.moveBy(-halfWidth, right);
        controlPoints[i*5+0] = topLeft;

        // bottom left
        Point3g bottomLeft = point;
        bottomLeft.moveBy(-halfHeight, normal);
        bottomLeft.moveBy(-halfWidth, right);
        controlPoints[i*5+1] = bottomLeft;

        // bottom right
        Point3g bottomRight = point;
        bottomRight.moveBy(-halfHeight, normal);
        bottomRight.moveBy(halfWidth, right);
        controlPoints[i*5+2] = bottomRight;

        // top right
        Point3g topRight = point;
        topRight.moveBy(halfHeight, normal);
        topRight.moveBy(halfWidth, right);
        controlPoints[i*5+3] = topRight;

        // top left (again)
        controlPoints[i*5+4] = topLeft;
    }

    QVector<GraphicsFloat> uKnots(trace.size() + degree + 1);
    for(int i = 0; i < uKnots.size(); i++){
        if(i < order)
            uKnots[i] = 0;
        else if(i > (uKnots.size() - 1 - order))
            uKnots[i] = uKnots.size()-2*order+1;
        else
            uKnots[i] = i-degree;
    }

    QVector<GraphicsFloat> vKnots(7);
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
