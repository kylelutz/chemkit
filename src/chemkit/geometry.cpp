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

// References:
//   - "Lecture Notes on Geometric Robustness" by Jonathan Richard Shewchuk
//   - "Measuring Space Filling Diagrams and Voids" by Herbert Edelsbrunner
//     and Ping Fu

#include "geometry.h"

#include "point3.h"
#include "vector.h"
#include "staticmatrix.h"
#include "staticvector.h"

namespace chemkit {

/// \ingroup chemkit
/// \brief The %chemkit::%geometry namespace contains various
///        construction and predicate functions for geometric
///        objects.
///
/// The functions in this namespace implement basic computational
/// geometry predicates and constructions. These are used as the
/// foundation for implementing higher-level geometric algorithms.
///
/// Geometric predicates are functions that return a true or false
/// result based on a geometric calculation. One example is the
/// planeOrientation() function which calculates whether a given
/// point lies above or below a plane.
///
/// Geometric constructions are functions that create a geometric
/// object and return information about it. One example is the
/// tetrahedronVolume() function which calculates the volume of a
/// tetrahedron given the positions of its four verticies.
namespace geometry {

// --- Constructions ------------------------------------------------------- //
Point3 circumcenter(const Point3 &a, const Point3 &b)
{
    Q_UNUSED(a);
    Q_UNUSED(b);

    return Point3();
}

Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c)
{
    Q_UNUSED(a);
    Q_UNUSED(b);
    Q_UNUSED(c);

    return Point3();
}

Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Q_UNUSED(a);
    Q_UNUSED(b);
    Q_UNUSED(c);
    Q_UNUSED(d);

    return Point3();
}

Float circumradius(const Point3 &a, const Point3 &b)
{
    Q_UNUSED(a);
    Q_UNUSED(b);

    return 0;
}

Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c)
{
    Q_UNUSED(a);
    Q_UNUSED(b);
    Q_UNUSED(c);

    return 0;
}

Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Q_UNUSED(a);
    Q_UNUSED(b);
    Q_UNUSED(c);
    Q_UNUSED(d);

    return 0;
}

/// Returns the position of a point orthogonal to the weighted
/// points (a, b).
Point3 orthocenter(const Point3 &a, const Point3 &b, Float wa, Float wb)
{
    double l = (1.0/2.0) - ((wa - wb) / (2.0 * (a - b).lengthSquared()));

    return a.scaled(l) + b.scaled(1.0 - l);
}

/// Returns the position of a point orthogonal to the weighted
/// points (a, b, c).
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc)
{
    // a1 = | a.y a.z 1 |
    //      | b.y b.z 1 |
    //      | c.y c.z 1 |
    double a1 = a.y()*b.z() + a.z()*c.y() + b.y()*c.z() -
                b.z()*c.y() - a.z()*b.y() - a.y()*c.z();

    // a2 = | a.z a.x 1 |
    //      | b.z b.x 1 |
    //      | c.z c.x 1 |
    double a2 = a.z()*b.x() + a.x()*c.z() + b.z()*c.x() -
                b.x()*c.z() - a.x()*b.z() - a.z()*c.x();

    // a3 = | a.x a.y 1 |
    //      | b.x b.y 1 |
    //      | c.x c.y 1 |
    double a3 = a.x()*b.y() + a.y()*c.x() + b.x()*c.y() -
                b.y()*c.x() - a.y()*b.x() - a.x()*c.y();

    // a4 = | a.x a.y a.z |
    //      | b.x b.y b.z |
    //      | c.x c.y c.z |
    double a4 = a.x()*b.y()*c.z() + a.y()*b.z()*c.x() + a.z()*b.x()*c.y() -
                a.z()*b.y()*c.x() - a.y()*b.x()*c.z() - a.x()*b.z()*c.y();

    double i0 = (1.0/2.0) * (a.dot(a) - wa);
    double j0 = (1.0/2.0) * (b.dot(b) - wb);
    double k0 = (1.0/2.0) * (c.dot(c) - wc);

    StaticMatrix<double, 4, 4> matrix;

    // d0 = | a.x a.y a.z 1 |
    //      | b.x b.y b.z 1 |
    //      | c.x c.y c.z 1 |
    //      | a1  a2  a3  0 |
    matrix << a.x(), a.y(), a.z(), 1,
              b.x(), b.y(), b.z(), 1,
              c.x(), c.y(), c.z(), 1,
              a1,    a2,    a3,    0;
    double d0 = matrix.determinant();

    // dx = | i0  a.y a.z 1 |
    //      | j0  b.y b.z 1 |
    //      | k0  c.y c.z 1 |
    //      | a1  a2  a3  0 |
    matrix << i0, a.y(), a.z(), 1,
              j0, b.y(), b.z(), 1,
              k0, c.y(), c.z(), 1,
              a4, a2,    a3,    0;
    double dx = matrix.determinant();

    // dy = | a.x i0  a.z 1 |
    //      | b.x j0  b.z 1 |
    //      | c.x k0  c.z 1 |
    //      | a1  a4  a3  0 |
    matrix << a.x(), i0, a.z(), 1,
              b.x(), j0, b.z(), 1,
              c.x(), k0, c.z(), 1,
              a1,    a4, a3,    0;
    double dy = matrix.determinant();

    // dz = | a.x a.y i0  1 |
    //      | b.x b.y j0  1 |
    //      | c.x c.y k0  1 |
    //      | a1  a2  a4  0 |
    matrix << a.x(), a.y(), i0, 1,
              b.x(), b.y(), j0, 1,
              c.x(), c.y(), k0, 1,
              a1,    a2,    a4, 0;
    double dz = matrix.determinant();

    return Point3(dx/d0, dy/d0, dz/d0);
}

/// Returns the position of a point orthogonal to the weighted
/// points (a, b, c, d).
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd)
{
    double i0 = (1.0/2.0) * (a.dot(a) - wa);
    double j0 = (1.0/2.0) * (b.dot(b) - wb);
    double k0 = (1.0/2.0) * (c.dot(c) - wc);
    double l0 = (1.0/2.0) * (d.dot(d) - wd);

    StaticMatrix<double, 4, 4> matrix;

    matrix << a.x(), a.y(), a.z(), 1,
              b.x(), b.y(), b.z(), 1,
              c.x(), c.y(), c.z(), 1,
              d.x(), d.y(), d.z(), 1;
    double d0 = matrix.determinant();

    matrix << i0, a.y(), a.z(), 1,
              j0, b.y(), b.z(), 1,
              k0, c.y(), c.z(), 1,
              l0, d.y(), d.z(), 1;
    double dx = matrix.determinant();

    matrix << a.x(), i0, a.z(), 1,
              b.x(), j0, b.z(), 1,
              c.x(), k0, c.z(), 1,
              d.x(), l0, d.z(), 1;
    double dy = matrix.determinant();

    matrix << a.x(), a.y(), i0, 1,
              b.x(), b.y(), j0, 1,
              c.x(), c.y(), k0, 1,
              d.x(), d.y(), l0, 1;
    double dz = matrix.determinant();

    return Point3(dx/d0, dy/d0, dz/d0);
}

/// Returns the radius of a point orthogonal to the edge with
/// weighted verticies (\p a, \p b).
Float orthoradius(const Point3 &a, const Point3 &b, Float wa, Float wb)
{
    GenericVector<double> ap = a - orthocenter(a, b, wa, wb);

    return ap.lengthSquared() - wa;
}

/// Returns the radius of a point orthogonal to the triangle
/// with weighted verticies (\p a, \p b, \p c).
Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc)
{
    GenericVector<double> r = a - c;
    GenericVector<double> s = b - c;

    double r2 = r.lengthSquared() + (wc - wa);
    double s2 = s.lengthSquared() + (wc - wb);

    // m1 = | (|r|^2 + (wc - wa)) r.y |
    //      | (|s|^2 + (wc - wb)) s.y |
    double m1 = r2 * s.y() - s2 * r.y();

    // m2 = | r.x (|r|^2 + (wc - wa)) |
    //      | s.x (|s|^2 + (wc - wb)) |
    double m2 = r.x() * s2 - s.x() * r2;

    // m3 = | (|r|^2 + (wc - wa)) r.z |
    //      | (|s|^2 + (wc - wb)) s.z |
    double m3 = r2 * s.z() - s2 * r.z();

    double A = r.cross(s).length() / 2.0;

    return (m1*m1 + m2*m2 + m3*m3) / (16.0 * A*A) - wc;
}

/// Returns the radius of a point orthogonal to the tetrahedron
/// with weighted verticies (\p a, \p b, \p c, \p d).
Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd)
{
    GenericVector<double> t = a - d;
    GenericVector<double> u = b - d;
    GenericVector<double> v = c - d;

    GenericVector<double> i = u.cross(v).scaled(t.lengthSquared() + (wd - wa));
    GenericVector<double> j = v.cross(t).scaled(u.lengthSquared() + (wd - wb));
    GenericVector<double> k = t.cross(u).scaled(v.lengthSquared() + (wd - wc));

    GenericVector<double> l = i + j + k;

    double V = tetrahedronVolume(a, b, c, d);

    return l.lengthSquared() / (144.0 * V*V) - wd;
}

/// Returns the area of the triangle with verticies (\p a, \p b,
/// \p c).
Float triangleArea(const Point3 &a, const Point3 &b, const Point3 &c)
{
    return (1.0/2.0) * (b - a).cross(c - a).length();
}

/// Returns the volume of the tetrahedron with verticies (\p a, \p b,
/// \p c, \p d).
Float tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    return (1.0/6.0) * planeOrientation(a, b, c, d);
}

// --- Predicates ---------------------------------------------------------- //
/// Returns a positive value if the point \p p is above the plane
/// defined by the points \p a, \p b, and \p c. Returns \c 0 if the
/// point lies on the plane. Returns a negative value if the point
/// lies below the plane.
///
/// Determination of the orientation is done by evaluating the sign
/// of the following determinant:
///
/** \f[
///   \left|
///   {
///     \begin{array}{cccc}
///       a_{x} & a_{y} & a_{z} & 1 \\
///       b_{x} & b_{y} & b_{z} & 1 \\
///       c_{x} & c_{y} & c_{z} & 1 \\
///       p_{x} & p_{y} & p_{z} & 1 \\
///     \end{array}
///   }
///   \right|
/// \f]
**/
Float planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p)
{
    GenericVector<double> t = a - p;
    GenericVector<double> u = b - p;
    GenericVector<double> v = c - p;

    return t.scalarTriple(u, v);
}

/// Returns a positive value if the point \p p is inside the sphere
/// made by points \p a, \p b, \p c, and \p d. Returns \c 0 if the
/// point is on the sphere. Returns a negative value if the point is
/// outside the sphere.
///
/// Determination of the sphere orientation is performed by
/// evaluating the sign of the following determinant:
///
/** \f[
///   \left|
///   {
///     \begin{array}{ccccc}
///       a_{x} & a_{y} & a_{z} & a_{x}^{2} + a_{y}^{2} + a_{z}^{2} & 1 \\
///       b_{x} & b_{y} & b_{z} & b_{x}^{2} + b_{y}^{2} + b_{z}^{2} & 1 \\
///       c_{x} & c_{y} & c_{z} & c_{x}^{2} + c_{y}^{2} + c_{z}^{2} & 1 \\
///       d_{x} & d_{y} & d_{z} & d_{x}^{2} + d_{y}^{2} + d_{z}^{2} & 1 \\
///       p_{x} & p_{y} & p_{z} & p_{x}^{2} + p_{y}^{2} + p_{z}^{2} & 1 \\
///     \end{array}
///   }
///   \right|
/// \f]
**/
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p)
{
    GenericVector<double> t = a - p;
    GenericVector<double> u = b - p;
    GenericVector<double> v = c - p;
    GenericVector<double> w = d - p;

    StaticMatrix<double, 4, 4> matrix;

    matrix << t.x(), t.y(), t.z(), t.dot(t),
              u.x(), u.y(), u.z(), u.dot(u),
              v.x(), v.y(), v.z(), v.dot(v),
              w.x(), w.y(), w.z(), w.dot(w);

    return matrix.determinant();
}

/// Returns a positive value if the weighted point \p p is inside
/// the sphere made by weighted points \p a, \p b, \p c, and \p d.
/// Returns \c 0 if the point is on the sphere. Returns a negative
/// value if the point is outside the sphere.
///
/// Determination of the weighted sphere orientation is performed
/// by evaluating the sign of the following determinant:
///
/** \f[
///   \left|
///   {
///     \begin{array}{ccccc}
///       a_{x} & a_{y} & a_{z} & a_{x}^{2} + a_{y}^{2} + a_{z}^{2} - wa & 1 \\
///       b_{x} & b_{y} & b_{z} & b_{x}^{2} + b_{y}^{2} + b_{z}^{2} - wb & 1 \\
///       c_{x} & c_{y} & c_{z} & c_{x}^{2} + c_{y}^{2} + c_{z}^{2} - wc & 1 \\
///       d_{x} & d_{y} & d_{z} & d_{x}^{2} + d_{y}^{2} + d_{z}^{2} - wd & 1 \\
///       p_{x} & p_{y} & p_{z} & p_{x}^{2} + p_{y}^{2} + p_{z}^{2} - wp & 1 \\
///     \end{array}
///   }
///   \right|
/// \f]
**/
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Float wa, Float wb, Float wc, Float wd, Float wp)
{
    GenericVector<double> t = a - p;
    GenericVector<double> u = b - p;
    GenericVector<double> v = c - p;
    GenericVector<double> w = d - p;

    StaticMatrix<double, 4, 4> matrix;

    matrix << t.x(), t.y(), t.z(), t.dot(t) - (wa - wp),
              u.x(), u.y(), u.z(), u.dot(u) - (wb - wp),
              v.x(), v.y(), v.z(), v.dot(v) - (wc - wp),
              w.x(), w.y(), w.z(), w.dot(w) - (wd - wp);

    return matrix.determinant();
}

} // end geometry namespace

} // end chemkit namespace
