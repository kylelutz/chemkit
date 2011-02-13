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

#include "alphashape.h"

#include "vector.h"
#include "geometry.h"
#include "staticmatrix.h"
#include "delaunaytriangulation.h"

namespace chemkit {

// === AlphaShapePrivate =================================================== //
class AlphaShapePrivate
{
    public:
        Float alphaValue;
        DelaunayTriangulation *triangulation;
};

// === AlphaShape ========================================================== //
/// \class AlphaShape alphashape.h chemkit/alphashape.h
/// \ingroup chemkit
/// \internal
/// \brief The AlphaShape class represents an alpha shape.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new alpha shape with \p points.
AlphaShape::AlphaShape(const QVector<Point> &points)
    : d(new AlphaShapePrivate)
{
    d->alphaValue = 0;
    d->triangulation = new DelaunayTriangulation(points);
}

/// Creates a new alpha shape with \p points and \p weights.
AlphaShape::AlphaShape(const QVector<Point> &points, const QVector<Float> &weights)
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
/// Returns the number of verticies in the alpha shape.
int AlphaShape::size() const
{
    return vertexCount();
}

/// Returns the position of \p vertex;
Point AlphaShape::position(int vertex) const
{
    return d->triangulation->position(vertex);
}

/// Returns the weight of \p vertex.
Float AlphaShape::weight(int vertex) const
{
    return d->triangulation->weight(vertex);
}

/// Sets the alpha value to \p alphaValue.
void AlphaShape::setAlphaValue(Float alphaValue)
{
    d->alphaValue = alphaValue;
}

/// Returns the alpha value.
Float AlphaShape::alphaValue() const
{
    return d->alphaValue;
}

// --- Simplicies ---------------------------------------------------------- //
/// Returns a list of verticies in the alpha shape.
QList<int> AlphaShape::verticies() const
{
    QList<int> verticies;

    return verticies;
}

/// Returns the number of verticies in the alpha shape.
int AlphaShape::vertexCount() const
{
    return verticies().size();
}

/// Returns a list of edges in the alpha shape.
QList<QVector<int> > AlphaShape::edges() const
{
    return d->triangulation->alphaShapeEdges(this);
}

/// Returns the number of edges in the alpha shape.
int AlphaShape::edgeCount() const
{
    return edges().size();
}

/// Returns a list of the triangles in the alpha shape.
QList<QVector<int> > AlphaShape::triangles() const
{
    return d->triangulation->alphaShapeTriangles(this);
}

/// Returns the number of triangles in the alpha shape.
int AlphaShape::triangleCount() const
{
    return triangles().size();
}

/// Returns a list of the tetrahedra in the alpha shape.
QList<QVector<int> > AlphaShape::tetrahedra() const
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
Float AlphaShape::volume() const
{
    Float volume = 0;

    foreach(const QVector<int> &tetrahedron, tetrahedra()){
        const Point &a = position(tetrahedron[0]);
        const Point &b = position(tetrahedron[1]);
        const Point &c = position(tetrahedron[2]);
        const Point &d = position(tetrahedron[3]);

        volume += chemkit::geometry::tetrahedronVolume(a, b, c, d);
    }

    return volume;
}

/// Returns the total surface area of the alpha shape.
Float AlphaShape::surfaceArea() const
{
    Float surfaceArea = 0;

//    foreach(const QVector<int> triangle, triangles(Regular | Singular)){
//        const Point &a = position(triangle[0]);
//        const Point &b = position(triangle[1]);
//        const Point &c = position(triangle[2]);
//
//        surfaceArea += chemkit::geometry::triangleArea(a, b, c);
//    }

    return surfaceArea;
}

Point AlphaShape::orthocenter(int i, int j) const
{
    const Point &a = position(i);
    const Point &b = position(j);

    Float wa = weight(i);
    Float wb = weight(j);

    return chemkit::geometry::orthocenter(a, b, wa, wb);
}

Point AlphaShape::orthocenter(int i, int j, int k) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);

    Float wa = weight(i);
    Float wb = weight(j);
    Float wc = weight(k);

    return chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
}

Point AlphaShape::orthocenter(int i, int j, int k, int l) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);
    const Point &d = position(l);

    Float wa = weight(i);
    Float wb = weight(j);
    Float wc = weight(k);
    Float wd = weight(l);

    return chemkit::geometry::orthocenter(a, b, c, d, wa, wb, wc, wd);
}

/// Returns the radius of the edge (a, b).
Float AlphaShape::orthoradius(int a, int b) const
{
    const Point &pa = position(a);
    const Point &pb = position(b);

    Float wa = weight(a);
    Float wb = weight(b);

    return chemkit::geometry::orthoradius(pa, pb, wa, wb);
}

/// Returns the radius of the triangle (a, b, c).
Float AlphaShape::orthoradius(int a, int b, int c) const
{
    const Point &pa = position(a);
    const Point &pb = position(b);
    const Point &pc = position(c);

    Float wa = weight(a);
    Float wb = weight(b);
    Float wc = weight(c);

    return chemkit::geometry::orthoradius(pa, pb, pc, wa, wb, wc);
}

Float AlphaShape::orthoradius(int a, int b, int c, int d) const
{
    const Point &pa = position(a);
    const Point &pb = position(b);
    const Point &pc = position(c);
    const Point &pd = position(d);

    Float wa = weight(a);
    Float wb = weight(b);
    Float wc = weight(c);
    Float wd = weight(d);

    return chemkit::geometry::orthoradius(pa, pb, pc, pd, wa, wb, wc, wd);
}

/// Returns \c true if the vertex i is attached to the vertex j.
bool AlphaShape::vertexAttached(int i, int j) const
{
    const Point &a = position(i);
    const Point &b = position(j);

    Float wa = weight(i);
    Float wb = weight(j);

    return (a - b).lengthSquared() + wa - wb < 0;
}

/// Returns \c true if the edge (i, j) is attached to vertex k.
bool AlphaShape::edgeAttached(int i, int j, int k) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);

    Float wa = weight(i);
    Float wb = weight(j);
    Float wc = weight(k);

    Point center = chemkit::geometry::orthocenter(a, b, wa, wb);
    Float radius = chemkit::geometry::orthoradius(a, b, wa, wb);

    return (center - c).lengthSquared() - radius - wc < 0;
}

/// Returns \c true if the triangle (\p i, \p j, \p k) is attached
/// to the vertex \p l.
bool AlphaShape::triangleAttached(int i, int j, int k, int l) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);
    const Point &d = position(l);

    Float wa = weight(i);
    Float wb = weight(j);
    Float wc = weight(k);
    Float wd = weight(l);

    Point center = chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
    Float radius = chemkit::geometry::orthoradius(a, b, c, wa, wb, wc);

    return (center - d).lengthSquared() - radius - wd < 0;
}

/// Returns \c true if the triangle (\p i, \p j, \p k) is attached
/// to either vertex \p l or vertex \p m.
bool AlphaShape::triangleAttached(int i, int j, int k, int l, int m) const
{
    const Point &a = position(i);
    const Point &b = position(j);
    const Point &c = position(k);
    const Point &d = position(l);
    const Point &e = position(m);

    Float wa = weight(i);
    Float wb = weight(j);
    Float wc = weight(k);
    Float wd = weight(l);
    Float we = weight(m);

    Point center = chemkit::geometry::orthocenter(a, b, c, wa, wb, wc);
    Float radius = chemkit::geometry::orthoradius(a, b, c, wa, wb, wc);

    if((center - d).lengthSquared() - radius - wd < 0){
        return true;
    }

    if((center - e).lengthSquared() - radius - we < 0){
        return true;
    }

    return false;
}

} // end chemkit namespace
