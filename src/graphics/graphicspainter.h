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

#ifndef CHEMKIT_GRAPHICSPAINTER_H
#define CHEMKIT_GRAPHICSPAINTER_H

#include "graphics.h"

namespace chemkit {

class Point3g;
class GraphicsVector;
class GraphicsMaterial;
class GraphicsVertexBuffer;
class GraphicsPainterPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsPainter
{
    public:
        // enumerations
        enum PrimitiveType {
            Triangles,
            TriangleStrip,
            TriangleFan,
            Lines,
            LineStrip,
            LineLoop,
            Points
        };

        // construction and destruction
        GraphicsPainter();
        ~GraphicsPainter();

        // drawing
        void draw(const GraphicsVertexBuffer *buffer, PrimitiveType type = Triangles);
        void drawSphere(GraphicsFloat radius);
        void drawSphere(const Point3g &center, GraphicsFloat radius);
        void drawCylinder(GraphicsFloat radius, GraphicsFloat length);
        void drawCylinder(const Point3g &a, const Point3g &b, GraphicsFloat radius);
        void drawCircle(GraphicsFloat radius);
        void drawCircle(const Point3g &center, GraphicsFloat radius, const GraphicsVector &normal);
        void drawTriangle(const Point3g &a, const Point3g &b, const Point3g &c);
        void drawRectangle(const Point3g &a, const Point3g &b, const Point3g &c, const Point3g &d);
        void drawSpline(const QList<Point3g> &points, GraphicsFloat radius, int order);
        void drawNurbsSurface(const QVector<Point3g> &controlPoints, const QVector<GraphicsFloat> &uKnots, const QVector<GraphicsFloat> &vKnots, int uOrder, int vOrder);
        void drawText(const QString &text, const QFont &font = QFont());
        void setColor(const QColor &color);
        void setMaterial(const GraphicsMaterial *material);

    private:
        GraphicsPainterPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPAINTER_H
