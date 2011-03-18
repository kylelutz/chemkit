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

#ifndef CHEMKIT_GRAPHICSVERTEXBUFFER_H
#define CHEMKIT_GRAPHICSVERTEXBUFFER_H

#include "graphics.h"

#include "point3g.h"
#include "graphicsvector.h"

namespace chemkit {

class GraphicsVertexBufferPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsVertexBuffer
{
    public:
        // enumerations
        enum Usage {
            Static,
            Dynamic,
            Stream
        };

        // construction and destruction
        GraphicsVertexBuffer();
        GraphicsVertexBuffer(const QVector<Point3g> &verticies);
        ~GraphicsVertexBuffer();

        // properties
        int size() const;
        bool isEmpty() const;
        void clear();

        // verticies
        void setVerticies(const QVector<Point3g> &verticies);
        QVector<Point3g> verticies() const;
        int vertexCount() const;

        // normals
        void setNormals(const QVector<GraphicsVector> &normals);
        QVector<GraphicsVector> normals() const;
        int normalCount() const;

        // indicies
        void setIndicies(const QVector<unsigned short> &indicies);
        QVector<unsigned short> indicies() const;
        int indexCount() const;

        // drawing
        void draw() const;
        void prepareToDraw() const;
        bool readyToDraw() const;

    private:
        GraphicsVertexBufferPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSVERTEXBUFFER_H
