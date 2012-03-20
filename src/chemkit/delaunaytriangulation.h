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

#ifndef CHEMKIT_DELAUNAYTRIANGULATION_H
#define CHEMKIT_DELAUNAYTRIANGULATION_H

#include "chemkit.h"

#include <vector>

#include <boost/array.hpp>

#include "point3.h"

namespace chemkit {

class AlphaShape;
class DelaunayTriangulationPrivate;

class CHEMKIT_EXPORT DelaunayTriangulation
{
public:
    // typedefs
    typedef boost::array<int, 2> Edge;
    typedef boost::array<int, 3> Triangle;

    // construction and destruction
    DelaunayTriangulation(const std::vector<Point3> &points);
    DelaunayTriangulation(const std::vector<Point3> &points, const std::vector<Real> &weights);
    ~DelaunayTriangulation();

    // properties
    int size() const;
    Point3 position(int vertex) const;
    Real weight(int vertex) const;
    bool isWeighted() const;

    // simplicies
    std::vector<int> vertices() const;
    int vertexCount() const;
    const std::vector<Edge>& edges() const;
    int edgeCount() const;
    const std::vector<Triangle>& triangles() const;
    int triangleCount() const;
    const std::vector<std::vector<int> >& tetrahedra() const;
    int tetrahedronCount() const;

    // geometry
    Real volume() const;
    Real surfaceArea() const;

private:
    // internal methods
    void triangulate(bool weighted);
    int location(const Point3 &point) const;
    void insertPoint(int index);
    std::vector<int> findContainingTetrahedra(int vertex) const;
    bool isExternal(int tetrahedron) const;

    // alpha shape
    const std::vector<Edge>& alphaShapeEdges(const AlphaShape *alphaShape) const;
    const std::vector<Triangle>& alphaShapeTriangles(const AlphaShape *alphaShape) const;
    const std::vector<std::vector<int> >& alphaShapeTetrahedra(const AlphaShape *alphaShape) const;
    void calculateAlphaShape(const AlphaShape *alphaShape) const;

    friend class AlphaShape;

private:
    DelaunayTriangulationPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_DELAUNAYTRIANGULATION_H
