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
CHEMKIT_EXPORT Float distance(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Float distanceSquared(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Float angle(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Float angleRadians(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Float torsionAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Float torsionAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Float wilsonAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Float wilsonAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Point3 midpoint(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Float circumradius(const Point3 &a, const Point3 &b);
CHEMKIT_EXPORT Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Float circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, Float wa, Float wb);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc);
CHEMKIT_EXPORT Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd);
CHEMKIT_EXPORT Float orthoradius(const Point3 &a, const Point3 &b, Float wa, Float wb);
CHEMKIT_EXPORT Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, Float wa, Float wb, Float wc);
CHEMKIT_EXPORT Float orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, Float wa, Float wb, Float wc, Float wd);
CHEMKIT_EXPORT Float triangleArea(const Point3 &a, const Point3 &b, const Point3 &c);
CHEMKIT_EXPORT Float tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
CHEMKIT_EXPORT Vector3 planeNormal(const Point3 &a, const Point3 &b, const Point3 &c);

// predicates
CHEMKIT_EXPORT Float planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p);
CHEMKIT_EXPORT Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p);
CHEMKIT_EXPORT Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Float wa, Float wb, Float wc, Float wd, Float wp);

} // end geometry namespace

} // end chemkit namespace

#endif // CHEMKIT_GEOMETRY_H
