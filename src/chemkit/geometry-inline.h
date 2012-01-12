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

#ifndef CHEMKIT_GEOMETRY_INLINE_H
#define CHEMKIT_GEOMETRY_INLINE_H

#include "geometry.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "constants.h"

namespace chemkit {

namespace geometry {

// --- Derivatives ------------------------------------------------------- //
/// Returns the gradient of the distance between points \p a and \p b.
inline boost::array<Vector3, 2> distanceGradient(const Point3 &a, const Point3 &b)
{
    boost::array<Vector3, 2> gradient;

    Real distance = chemkit::geometry::distance(a, b);

    gradient[0] = (a - b) / distance;
    gradient[1] = -gradient[0];

    return gradient;
}

/// Returns the gradient of the angle between points \p a, \p b
/// and \p c.
inline boost::array<Vector3, 3> angleGradient(const Point3 &a, const Point3 &b, const Point3 &c)
{
    boost::array<Vector3, 3> gradient = chemkit::geometry::angleGradientRadians(a, b, c);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the angle between points \p a, \p b
/// and \p c.
inline boost::array<Vector3, 3> angleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c)
{
    boost::array<Vector3, 3> gradient;

    Real theta = chemkit::geometry::angleRadians(a, b, c);

    Real rab = chemkit::geometry::distance(a, b);
    Real rbc = chemkit::geometry::distance(b, c);

    gradient[0] = ((((c - b) * rab) - (a - b) * ((b - a).dot(b - c) / rab)) / (pow(rab, 2) * rbc)) / -sin(theta);
    gradient[1] = ((((b - c) + (b - a)) * (rab * rbc) - (((b - a) * (rbc/rab) + (b - c) * (rab/rbc)) * (b - a).dot(b - c))) / pow(rab * rbc, 2)) / -sin(theta);
    gradient[2] = -gradient[0] - gradient[1];

    return gradient;
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> torsionAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient = chemkit::geometry::torsionAngleGradientRadians(a, b, c, d);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> torsionAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient;

    Real phi = chemkit::geometry::torsionAngleRadians(a, b, c, d);

    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 bd = d - b;
    Vector3 cb = b - c;
    Vector3 cd = d - c;

    Vector3 m = ab.cross(cb);
    Vector3 n = cb.cross(cd);

    Vector3 p = ((n / (m.norm() * n.norm())) - ((m / m.squaredNorm()) * cos(phi)));
    Vector3 q = ((m / (m.norm() * n.norm())) - ((n / n.squaredNorm()) * cos(phi)));

    gradient[0] = cb.cross(p) * (1.0 / sin(phi));
    gradient[1] = (ac.cross(p) - cd.cross(q)) * (1.0 / sin(phi));
    gradient[2] = (bd.cross(q) - ab.cross(p)) * (1.0 / sin(phi));
    gradient[3] = cb.cross(q) * (1.0 / sin(phi));

    return gradient;
}

/// Returns the gradient of the wilson angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> wilsonAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient = chemkit::geometry::wilsonAngleGradientRadians(a, b, c, d);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the wilson angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> wilsonAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Vector3 ba = a - b;
    Vector3 bc = c - b;
    Vector3 bd = d - b;

    Real rba = ba.norm();
    Real rbc = bc.norm();
    Real rbd = bd.norm();

    ba.normalize();
    bc.normalize();
    bd.normalize();

    Real theta = acos(ba.dot(bc));

    Real w = chemkit::geometry::wilsonAngleRadians(a, b, c, d);

    boost::array<Vector3, 4> gradient;

    gradient[0] = ((bd.cross(bc) / (cos(w) * sin(theta)) - (ba - bc * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rba;
    gradient[2] = ((ba.cross(bd) / (cos(w) * sin(theta)) - (bc - ba * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rbc;
    gradient[3] = (bc.cross(ba) / (cos(w) * sin(theta)) - bd * tan(w)) / rbd;
    gradient[1] = -(gradient[0] + gradient[2] + gradient[3]);

    return gradient;
}

// --- Transforms ---------------------------------------------------------- //
template<typename T>
inline Eigen::Matrix<T, 3, 1> rotate(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle)
{
    return rotateRadians<T>(vector, axis, angle * chemkit::constants::DegreesToRadians);
}

template<typename T>
inline Eigen::Matrix<T, 3, 1> rotateRadians(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle)
{
    return Eigen::Transform<T, 3, 3>(Eigen::AngleAxis<T>(angle, axis)) * vector;
}

} // end geometry namespace

} // end chemkit namespace

#endif // CHEMKIT_GEOMETRY_INLINE_H
