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

#ifndef CHEMKIT_GEOMETRY_H
#define CHEMKIT_GEOMETRY_H

#include "chemkit.h"

#include "point3.h"

namespace chemkit {

namespace geometry {

// constructions
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

// predicates
Float planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p);
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p);
Float sphereOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, const Point3 &p, Float wa, Float wb, Float wc, Float wd, Float wp);

} // end geometry namespace

} // end chemkit namespace

#endif // CHEMKIT_GEOMETRY_H
