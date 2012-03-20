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

#include "alphashape.h"

#include "foreach.h"
#include "vector3.h"
#include "geometry.h"

namespace chemkit {

// === AlphaShapePrivate =================================================== //
class AlphaShapePrivate
{
public:
    Real alphaValue;
    DelaunayTriangulation *triangulation;
};

// === AlphaShape ========================================================== //
/// \class AlphaShape alphashape.h chemkit/alphashape.h
/// \ingroup chemkit
/// \internal
/// \brief The AlphaShape class represents an alpha shape.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new alpha shape with \p points.
AlphaShape::AlphaShape(const std::vector<Point3> &points)
    : d(new AlphaShapePrivate)
{
    d->alphaValue = 0;
    d->triangulation = new DelaunayTriangulation(points);
}

/// Creates a new alpha shape with \p points and \p weights.
AlphaShape::AlphaShape(const std::vector<Point3> &points, const std::vector<Real> &weights)
    : d(new AlphaShapePrivate)
{
    d->alphaValue = 0;
    d->triangulation = new DelaunayTriangulation(points, weights);
}

/// Destroys the alpha shape object;
AlphaShape::~AlphaShape()
{
    delete d->triangulation;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of vertices in the alpha shape.
int AlphaShape::size() const
{
    return vertexCount();
}

/// Returns the position of \p vertex;
Point3 AlphaShape::position(int vertex) const
{
    return d->triangulation->position(vertex);
}

/// Returns the weight of \p vertex.
Real AlphaShape::weight(int vertex) const
{
    return d->triangulation->weight(vertex);
}

/// Sets the alpha value to \p alphaValue.
void AlphaShape::setAlphaValue(Real alphaValue)
{
    d->alphaValue = alphaValue;
}

/// Returns the alpha value.
Real AlphaShape::alphaValue() const
{
    return d->alphaValue;
}

// --- Simplicies ---------------------------------------------------------- //
/// Returns a list of vertices in the alpha shape.
std::vector<int> AlphaShape::vertices() const
{
    std::vector<int> vertices;

    return vertices;
}

/// Returns the number of vertices in the alpha shape.
int AlphaShape::vertexCount() const
{
    return vertices().size();
}

/// Returns a list of edges in the alpha shape.
const std::vector<AlphaShape::Edge>& AlphaShape::edges() const
{
    return d->triangulation->alphaShapeEdges(this);
}

/// Returns the number of edges in the alpha shape.
int AlphaShape::edgeCount() const
{
    return edges().size();
}

/// Returns a list of the triangles in the alpha shape.
const std::vector<AlphaShape::Triangle>& AlphaShape::triangles() const
{
    return d->triangulation->alphaShapeTriangles(this);
}

/// Returns the number of triangles in the alpha shape.
int AlphaShape::triangleCount() const
{
    return triangles().size();
}

/// Returns a list of the tetrahedra in the alpha shape.
const std::vector<std::vector<int> >& AlphaShape::tetrahedra() const
{
    return d->triangulation->alphaShapeTetrahedra(this);
}

/// Returns the number of tetrahedra in the alpha shape.
int AlphaShape::tetrahedronCount() const
{
    return tetrahedra().size();
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the total volume of the alpha shape.
Real AlphaShape::volume() const
{
    Real volume = 0;

    foreach(const std::vector<int> &tetrahedron, tetrahedra()){
        const Point3 &a = position(tetrahedron[0]);
        const Point3 &b = position(tetrahedron[1]);
        const Point3 &c = position(tetrahedron[2]);
        const Point3 &d = position(tetrahedron[3]);

        volume += chemkit::geometry::tetrahedronVolume(a, b, c, d);
    }

    return volume;
}

/// Returns the total surface area of the alpha shape.
Real AlphaShape::surfaceArea() const
{
    Real surfaceArea = 0;

//    foreach(const std::vector<int> triangle, triangles(Regular | Singular)){
//        const Point &a = position(triangle[0]);
//        const Point &b = position(triangle[1]);
//        const Point &c = position(triangle[2]);
//
//        surfaceArea += chemkit::geometry::triangleArea(a, b, c);
//    }

    return surfaceArea;
}

Point3 AlphaShape::orthocenter(int i, int j) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);

    Real wa = weight(i);
    Real wb = weight(j);

    return chemkit::geometry::orthocenter(a, b, wa, wb);
}

Point3 AlphaShape::orthocenter(int i, int j, int k) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);

    Real wa = weight(i);
    Real wb = weight(j);
    Real wc = weight(k);

    return chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
}

Point3 AlphaShape::orthocenter(int i, int j, int k, int l) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);
    const Point3 &d = position(l);

    Real wa = weight(i);
    Real wb = weight(j);
    Real wc = weight(k);
    Real wd = weight(l);

    return chemkit::geometry::orthocenter(a, b, c, d, wa, wb, wc, wd);
}

/// Returns the radius of the edge (a, b).
Real AlphaShape::orthoradius(int a, int b) const
{
    const Point3 &pa = position(a);
    const Point3 &pb = position(b);

    Real wa = weight(a);
    Real wb = weight(b);

    return chemkit::geometry::orthoradius(pa, pb, wa, wb);
}

/// Returns the radius of the triangle (a, b, c).
Real AlphaShape::orthoradius(int a, int b, int c) const
{
    const Point3 &pa = position(a);
    const Point3 &pb = position(b);
    const Point3 &pc = position(c);

    Real wa = weight(a);
    Real wb = weight(b);
    Real wc = weight(c);

    return chemkit::geometry::orthoradius(pa, pb, pc, wa, wb, wc);
}

Real AlphaShape::orthoradius(int a, int b, int c, int d) const
{
    const Point3 &pa = position(a);
    const Point3 &pb = position(b);
    const Point3 &pc = position(c);
    const Point3 &pd = position(d);

    Real wa = weight(a);
    Real wb = weight(b);
    Real wc = weight(c);
    Real wd = weight(d);

    return chemkit::geometry::orthoradius(pa, pb, pc, pd, wa, wb, wc, wd);
}

/// Returns \c true if the vertex i is attached to the vertex j.
bool AlphaShape::vertexAttached(int i, int j) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);

    Real wa = weight(i);
    Real wb = weight(j);

    return (a - b).squaredNorm() + wa - wb < 0;
}

/// Returns \c true if the edge (i, j) is attached to vertex k.
bool AlphaShape::edgeAttached(int i, int j, int k) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);

    Real wa = weight(i);
    Real wb = weight(j);
    Real wc = weight(k);

    Point3 center = chemkit::geometry::orthocenter(a, b, wa, wb);
    Real radius = chemkit::geometry::orthoradius(a, b, wa, wb);

    return (center - c).squaredNorm() - radius - wc < 0;
}

/// Returns \c true if the triangle (\p i, \p j, \p k) is attached
/// to the vertex \p l.
bool AlphaShape::triangleAttached(int i, int j, int k, int l) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);
    const Point3 &d = position(l);

    Real wa = weight(i);
    Real wb = weight(j);
    Real wc = weight(k);
    Real wd = weight(l);

    Point3 center = chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
    Real radius = chemkit::geometry::orthoradius(a, b, c, wa, wb, wc);

    return (center - d).squaredNorm() - radius - wd < 0;
}

/// Returns \c true if the triangle (\p i, \p j, \p k) is attached
/// to either vertex \p l or vertex \p m.
bool AlphaShape::triangleAttached(int i, int j, int k, int l, int m) const
{
    const Point3 &a = position(i);
    const Point3 &b = position(j);
    const Point3 &c = position(k);
    const Point3 &d = position(l);
    const Point3 &e = position(m);

    Real wa = weight(i);
    Real wb = weight(j);
    Real wc = weight(k);
    Real wd = weight(l);
    Real we = weight(m);

    Point3 center = chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
    Real radius = chemkit::geometry::orthoradius(a, b, c, wa, wb, wc);

    if((center - d).squaredNorm() - radius - wd < 0){
        return true;
    }

    if((center - e).squaredNorm() - radius - we < 0){
        return true;
    }

    return false;
}

} // end chemkit namespace
