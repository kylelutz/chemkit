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

#include <vector>

#include "point3.h"

namespace chemkit {

class AlphaShape;
class DelaunayTriangulationPrivate;

class CHEMKIT_EXPORT DelaunayTriangulation
{
    public:
        // construction and destruction
        DelaunayTriangulation(const std::vector<Point3> &points);
        DelaunayTriangulation(const std::vector<Point3> &points, const std::vector<Float> &weights);
        ~DelaunayTriangulation();

        // properties
        int size() const;
        Point3 position(int vertex) const;
        Float weight(int vertex) const;
        bool isWeighted() const;

        // simplicies
        std::vector<int> verticies() const;
        int vertexCount() const;
        const std::vector<std::vector<int> >& edges() const;
        int edgeCount() const;
        const std::vector<std::vector<int> >& triangles() const;
        int triangleCount() const;
        const std::vector<std::vector<int> >& tetrahedra() const;
        int tetrahedronCount() const;

        // geometry
        Float volume() const;
        Float surfaceArea() const;

    private:
        // internal methods
        void triangulate(bool weighted);
        int location(const Point3 &point) const;
        void insertPoint(int index);
        std::vector<int> findContainingTetrahedra(int vertex) const;
        bool isExternal(int tetrahedron) const;

        // alpha shape
        const std::vector<std::vector<int> >& alphaShapeEdges(const AlphaShape *alphaShape) const;
        const std::vector<std::vector<int> >& alphaShapeTriangles(const AlphaShape *alphaShape) const;
        const std::vector<std::vector<int> >& alphaShapeTetrahedra(const AlphaShape *alphaShape) const;
        void calculateAlphaShape(const AlphaShape *alphaShape) const;

        friend class AlphaShape;

    private:
        DelaunayTriangulationPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_DELAUNAYTRIANGULATION_H
