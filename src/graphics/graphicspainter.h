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

#ifndef CHEMKIT_GRAPHICSPAINTER_H
#define CHEMKIT_GRAPHICSPAINTER_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>

namespace chemkit {

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
        void drawSphere(float radius);
        void drawSphere(const Point3f &center, float radius);
        void drawCylinder(float radius, float length);
        void drawCylinder(const Point3f &a, const Point3f &b, float radius);
        void drawCircle(float radius);
        void drawCircle(const Point3f &center, float radius, const Vector3f &normal);
        void drawTriangle(const Point3f &a, const Point3f &b, const Point3f &c);
        void drawRectangle(const Point3f &a, const Point3f &b, const Point3f &c, const Point3f &d);
        void drawSpline(const QList<Point3f> &points, float radius, int order);
        void drawNurbsSurface(const QVector<Point3f> &controlPoints, const QVector<float> &uKnots, const QVector<float> &vKnots, int uOrder, int vOrder);
        void drawText(const QString &text, const QFont &font = QFont());
        void setColor(const QColor &color);
        void setMaterial(const GraphicsMaterial *material);

    private:
        GraphicsPainterPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSPAINTER_H
