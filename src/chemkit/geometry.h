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

#ifndef CHEMKIT_GEOMETRY_H
#define CHEMKIT_GEOMETRY_H

#include "chemkit.h"

#include <boost/array.hpp>

#include "point3.h"
#include "vector3.h"

namespace chemkit {

namespace geometry {

// constructions
CHEMKIT_EXPORT Real distance(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Real distanceSquared(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Real angle(const Vector3 &a, const Vector3 &b);
CHEMKIT_EXPORT Real angleRadians(const Vector3 &a, const Vector3 &b);
CHEMKIT_EXPORT Real angle(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Real angleRadians(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Real torsionAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Real torsionAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Real wilsonAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Real wilsonAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Point3 midpoint(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Real circumradius(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Real circumradius(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Real circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, Real wa, Real wb);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, Real wa, Real wb, Real wc);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Real wa, Real wb, Real wc, Real wd);
CHEMKIT_EXPORT Real orthoradius(const Point3 &a, const Point3 &b, Real wa, Real wb);
CHEMKIT_EXPORT Real orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, Real wa, Real wb, Real wc);
CHEMKIT_EXPORT Real orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Real wa, Real wb, Real wc, Real wd);
CHEMKIT_EXPORT Real triangleArea(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Real tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Vector3 planeNormal(const Point3 &a, const Point3 &b, const Point3 &c);

// derivatives
inline boost::array<Vector3, 2> distanceGradient(const Point3 &a, const Point3 &b);
inline boost::array<Vector3, 3> angleGradient(const Point3 &a, const Point3 &b, const Point3 &c);
inline boost::array<Vector3, 3> angleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c);
inline boost::array<Vector3, 4> torsionAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> torsionAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> wilsonAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> wilsonAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);

// predicates
CHEMKIT_EXPORT Real planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p);
CHEMKIT_EXPORT Real sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p);
CHEMKIT_EXPORT Real sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Real wa, Real wb, Real wc, Real wd, Real wp);

// transforms
template<typename T> Eigen::Matrix<T, 3, 1> rotate(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle);
template<typename T> Eigen::Matrix<T, 3, 1> rotateRadians(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle);

} // end geometry namespace

} // end chemkit namespace

#include "geometry-inline.h"

#endif // CHEMKIT_GEOMETRY_H
