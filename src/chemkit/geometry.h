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

#include "point3.h"
#include "vector3.h"

namespace chemkit {

namespace geometry {

// constructions
Float distance(const Point3 &a, const Point3 &b);
Float distanceSquared(const Point3 &a, const Point3 &b);
Float angle(const Point3 &a, const Point3 &b, const Point3 &c);
Float angleRadians(const Point3 &a, const Point3 &b, const Point3 &c);
Float torsionAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Float torsionAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Float wilsonAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Float wilsonAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Point3 midpoint(const Point3 &a, const Point3 &b);
Point3 circumcenter(const Point3 &a, const Point3 &b);
Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c);
Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Float circumradius(const Point3 &a, const Point3 &b);
Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c);
Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Point3 orthocenter(const Point3 &a, const Point3 &b, Float wa, Float wb);
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc);
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd);
Float orthoradius(const Point3 &a, const Point3 &b, Float wa, Float wb);
Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc);
Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd);
Float triangleArea(const Point3 &a, const Point3 &b, const Point3 &c);
Float tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Vector3 planeNormal(const Point3 &a, const Point3 &b, const Point3 &c);

// predicates
Float planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p);
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p);
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Float wa, Float wb, Float wc, Float wd, Float wp);

} // end geometry namespace

} // end chemkit namespace

#endif // CHEMKIT_GEOMETRY_H
