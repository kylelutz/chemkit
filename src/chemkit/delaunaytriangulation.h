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

#ifndef CHEMKIT_DELAUNAYTRIANGULATION_H
#define CHEMKIT_DELAUNAYTRIANGULATION_H

#include "chemkit.h"

#include <QtCore>

#include "point.h"

namespace chemkit {

class AlphaShape;
class DelaunayTriangulationPrivate;

class CHEMKIT_EXPORT DelaunayTriangulation
{
    public:
        // construction and destruction
        DelaunayTriangulation(const QVector<Point> &points);
        DelaunayTriangulation(const QVector<Point> &points, const QVector<Float> &weights);
        ~DelaunayTriangulation();

        // properties
        int size() const;
        Point position(int vertex) const;
        Float weight(int vertex) const;
        bool isWeighted() const;

        // simplicies
        QList<int> verticies() const;
        int vertexCount() const;
        QList<QVector<int> > edges() const;
        int edgeCount() const;
        QList<QVector<int> > triangles() const;
        int triangleCount() const;
        QList<QVector<int> > tetrahedra() const;
        int tetrahedronCount() const;

        // geometry
        Float volume() const;
        Float surfaceArea() const;

    private:
        // internal methods
        void triangulate(bool weighted);
        int location(const Point &point) const;
        void insertPoint(int index);
        QList<int> findContainingTetrahedra(int vertex) const;
        bool isExternal(int tetrahedron) const;

        // alpha shape
        QList<QVector<int> > alphaShapeEdges(const AlphaShape *alphaShape) const;
        QList<QVector<int> > alphaShapeTriangles(const AlphaShape *alphaShape) const;
        QList<QVector<int> > alphaShapeTetrahedra(const AlphaShape *alphaShape) const;
        void calculateAlphaShape(const AlphaShape *alphaShape) const;

        friend class AlphaShape;

    private:
        DelaunayTriangulationPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_DELAUNAYTRIANGULATION_H
