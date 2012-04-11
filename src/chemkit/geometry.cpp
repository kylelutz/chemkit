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

// References:
//   - "Lecture Notes on Geometric Robustness" by Jonathan Richard Shewchuk
//   - "Measuring Space Filling Diagrams and Voids" by Herbert Edelsbrunner
//     and Ping Fu

#include "geometry.h"

#include "point3.h"
#include "vector3.h"

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
/// tetrahedron given the positions of its four vertices.
namespace geometry {

// --- Constructions ------------------------------------------------------- //
/// Returns the distance between the points \p a and \p b.
Real distance(const Point3 &a, const Point3 &b)
{
    return (a - b).norm();
}

/// Returns the square of the distance between the points \p a and
/// \p b.
Real distanceSquared(const Point3 &a, const Point3 &b)
{
    return (a - b).squaredNorm();
}

/// Returns the angle between the vectors \p a and \p b. Angle is
/// in Degrees.
Real angle(const Vector3 &a, const Vector3 &b)
{
    return angleRadians(a, b) * chemkit::constants::RadiansToDegrees;
}

/// Returns the angle between the vectors \p a and \p b. Angle is
/// in Radians.
Real angleRadians(const Vector3 &a, const Vector3 &b)
{
    return acos(a.dot(b) / (a.norm() * b.norm()));
}

/// Returns the angle between the vectors (\p a, \p b) and (\p b,
/// \p c). Angle is in Degrees.
Real angle(const Point3 &a, const Point3 &b, const Point3 &c)
{
    return angleRadians(a, b, c) * chemkit::constants::RadiansToDegrees;
}

/// Returns the ang le between the vectors (\p a, \p b) and (\p b,
/// \p c). Angle is Radians.
Real angleRadians(const Point3 &a, const Point3 &b, const Point3 &c)
{
    Point3 ab = b - a;
    Point3 cb = b - c;

    return acos(ab.dot(cb) / (ab.norm() * cb.norm()));
}

/// Returns the torsion angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
Real torsionAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    return torsionAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

/// Returns the torsion angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
Real torsionAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Point3 ab = b - a;
    Point3 bc = c - b;
    Point3 cd = d - c;

    return atan2(bc.norm() * ab.dot(bc.cross(cd)), ab.cross(bc).dot(bc.cross(cd)));
}

/// Returns the wilson angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Degrees.
Real wilsonAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    return wilsonAngleRadians(a, b, c, d) * chemkit::constants::RadiansToDegrees;
}

/// Returns the wilson angle between points \p a, \p b, \p c, and
/// \p d. Angle is in Radians.
Real wilsonAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Vector3 normal = ((b - a).cross(c - b)).normalized();

    Point3 bd = d - b;
    bd.normalize();

    Real angle = acos(bd.dot(normal));

    return (chemkit::constants::Pi * 0.5) - angle;
}

/// Returns the midpoint between \p a and \p b.
Point3 midpoint(const Point3 &a, const Point3 &b)
{
    return (a + b) * 0.5;
}

Point3 circumcenter(const Point3 &a, const Point3 &b)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);

    return Point3();
}

Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);

    return Point3();
}

Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);
    CHEMKIT_UNUSED(d);

    return Point3();
}

Real circumradius(const Point3 &a, const Point3 &b)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);

    return 0;
}

Real circumradius(const Point3 &a, const Point3 &b, const Point3 &c)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);

    return 0;
}

Real circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);
    CHEMKIT_UNUSED(d);

    return 0;
}

/// Returns the position of a point orthogonal to the weighted
/// points (a, b).
Point3 orthocenter(const Point3 &a, const Point3 &b, Real wa, Real wb)
{
    double l = (1.0/2.0) - ((wa - wb) / (2.0 * (a - b).squaredNorm()));

    return (a * l) + (b * (1 - l));
}

/// Returns the position of a point orthogonal to the weighted
/// points (a, b, c).
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, Real wa, Real wb, Real wc)
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

    Eigen::Matrix<double, 4, 4> matrix;

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
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Real wa, Real wb, Real wc, Real wd)
{
    double i0 = (1.0/2.0) * (a.dot(a) - wa);
    double j0 = (1.0/2.0) * (b.dot(b) - wb);
    double k0 = (1.0/2.0) * (c.dot(c) - wc);
    double l0 = (1.0/2.0) * (d.dot(d) - wd);

    Eigen::Matrix<double, 4, 4> matrix;

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
/// weighted vertices (\p a, \p b).
Real orthoradius(const Point3 &a, const Point3 &b, Real wa, Real wb)
{
    Vector3d ap = a.cast<double>() - orthocenter(a, b, wa, wb).cast<double>();

    return ap.squaredNorm() - wa;
}

/// Returns the radius of a point orthogonal to the triangle
/// with weighted vertices (\p a, \p b, \p c).
Real orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, Real wa, Real wb, Real wc)
{
    Vector3d r = a.cast<double>() - c.cast<double>();
    Vector3d s = b.cast<double>() - c.cast<double>();

    double r2 = r.squaredNorm() + (wc - wa);
    double s2 = s.squaredNorm() + (wc - wb);

    // m1 = | (|r|^2 + (wc - wa)) r.y |
    //      | (|s|^2 + (wc - wb)) s.y |
    double m1 = r2 * s.y() - s2 * r.y();

    // m2 = | r.x (|r|^2 + (wc - wa)) |
    //      | s.x (|s|^2 + (wc - wb)) |
    double m2 = r.x() * s2 - s.x() * r2;

    // m3 = | (|r|^2 + (wc - wa)) r.z |
    //      | (|s|^2 + (wc - wb)) s.z |
    double m3 = r2 * s.z() - s2 * r.z();

    double A = r.cross(s).norm() / 2.0;

    return (m1*m1 + m2*m2 + m3*m3) / (16.0 * A*A) - wc;
}

/// Returns the radius of a point orthogonal to the tetrahedron
/// with weighted vertices (\p a, \p b, \p c, \p d).
Real orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Real wa, Real wb, Real wc, Real wd)
{
    Vector3d t = a.cast<double>() - d.cast<double>();
    Vector3d u = b.cast<double>() - d.cast<double>();
    Vector3d v = c.cast<double>() - d.cast<double>();

    Vector3d i = u.cross(v) * (t.squaredNorm() + (wd - wa));
    Vector3d j = v.cross(t) * (u.squaredNorm() + (wd - wb));
    Vector3d k = t.cross(u) * (v.squaredNorm() + (wd - wc));

    Vector3d l = i + j + k;

    double V = tetrahedronVolume(a, b, c, d);

    return l.squaredNorm() / (144.0 * V*V) - wd;
}

/// Returns the area of the triangle with vertices (\p a, \p b,
/// \p c).
Real triangleArea(const Point3 &a, const Point3 &b, const Point3 &c)
{
    return (1.0/2.0) * (b - a).cross(c - a).norm();
}

/// Returns the volume of the tetrahedron with vertices (\p a, \p b,
/// \p c, \p d).
Real tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    return (1.0/6.0) * planeOrientation(a, b, c, d);
}

/// Returns the normal vector to the plane made by the points \p a,
/// \p b, \p c.
Vector3 planeNormal(const Point3 &a, const Point3 &b, const Point3 &c)
{
    return ((b - a).cross(c - b)).normalized();
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
Real planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p)
{
    Vector3d t = a.cast<double>() - p.cast<double>();
    Vector3d u = b.cast<double>() - p.cast<double>();
    Vector3d v = c.cast<double>() - p.cast<double>();

    return t.dot(u.cross(v));
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
Real sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p)
{
    Vector3d t = a.cast<double>() - p.cast<double>();
    Vector3d u = b.cast<double>() - p.cast<double>();
    Vector3d v = c.cast<double>() - p.cast<double>();
    Vector3d w = d.cast<double>() - p.cast<double>();

    Eigen::Matrix<double, 4, 4> matrix;

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
Real sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Real wa, Real wb, Real wc, Real wd, Real wp)
{
    Vector3d t = a.cast<double>() - p.cast<double>();
    Vector3d u = b.cast<double>() - p.cast<double>();
    Vector3d v = c.cast<double>() - p.cast<double>();
    Vector3d w = d.cast<double>() - p.cast<double>();

    Eigen::Matrix<double, 4, 4> matrix;

    matrix << t.x(), t.y(), t.z(), t.dot(t) - (wa - wp),
              u.x(), u.y(), u.z(), u.dot(u) - (wb - wp),
              v.x(), v.y(), v.z(), v.dot(v) - (wc - wp),
              w.x(), w.y(), w.z(), w.dot(w) - (wd - wp);

    return matrix.determinant();
}

} // end geometry namespace

} // end chemkit namespace
