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

#ifndef CHEMKIT_FORCEFIELDCALCULATION_INLINE_H
#define CHEMKIT_FORCEFIELDCALCULATION_INLINE_H

#include "forcefieldcalculation.h"

#include "point.h"
#include "forcefieldatom.h"

namespace chemkit {

/// Returns the distance between atoms \p a and \p b. Distance is in
/// Angstroms.
inline Float ForceFieldCalculation::distance(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return a->position().distance(b->position());
}

/// Returns the gradient of the distance between atoms \p a and \p b.
inline QVector<Vector> ForceFieldCalculation::distanceGradient(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return distanceGradient(a->position(), b->position());
}

/// Returns the gradient of the distance between points \p a and \p b.
inline QVector<Vector> ForceFieldCalculation::distanceGradient(const Point &a, const Point &b) const
{
    QVector<Vector> gradient(2);

    gradient[0] = (a - b) / a.distance(b);
    gradient[1] = -gradient[0];

    return gradient;
}

/// Returns the bond angle between atoms \p a, \p b and \p c. The
/// angle is in degrees.
inline Float ForceFieldCalculation::bondAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return Point::angle(a->position(), b->position(), c->position());
}

/// Returns the bond angle between atoms \p a, \p b and \p c. The
/// angle is in radians.
inline Float ForceFieldCalculation::bondAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return Point::angleRadians(a->position(), b->position(), c->position());
}

/// Returns the gradient of the bond angle between atoms \p a, \p b
/// and \p c.
inline QVector<Vector> ForceFieldCalculation::bondAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    QVector<Vector> gradient = bondAngleGradientRadians(a, b, c);

    for(int i = 0; i < gradient.size(); i++){
        gradient[i].scale(chemkit::constants::RadiansToDegrees);
    }

    return gradient;
}

/// Returns the gradient of the bond angle between atoms \p a, \p b
/// and \p c.
inline QVector<Vector> ForceFieldCalculation::bondAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return bondAngleGradientRadians(a->position(), b->position(), c->position());
}

/// Returns the gradient of the bond angle between points \p a, \p b
/// and \p c.
inline QVector<Vector> ForceFieldCalculation::bondAngleGradientRadians(const Point &a, const Point &b, const Point &c) const
{
    QVector<Vector> gradient(3);

    Float theta = Point::angleRadians(a, b, c);

    Float rab = a.distance(b);
    Float rbc = b.distance(c);

    gradient[0] = ((((c - b) * rab) - (a - b) * ((b - a).dot(b - c) / rab)) / (pow(rab, 2) * rbc)) / -sin(theta);
    gradient[1] = ((((b - c) + (b - a)) * (rab * rbc) - (((b - a) * (rbc/rab) + (b - c) * (rab/rbc)) * (b - a).dot(b - c))) / pow(rab * rbc, 2)) / -sin(theta);
    gradient[2] = -gradient[0] - gradient[1];

    return gradient;
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between the atoms \p a, \p b, \p c, and \p d. The angle is in
/// degrees.
inline Float ForceFieldCalculation::torsionAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point::torsionAngle(a->position(), b->position(), c->position(), d->position());
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between the atoms \p a, \p b, \p c, and \p d. The angle is in
/// radians.
inline Float ForceFieldCalculation::torsionAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point::torsionAngleRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the torsion angle between the atoms \p a,
/// \p b, \p c, and \p d.
inline QVector<Vector> ForceFieldCalculation::torsionAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    QVector<Vector> gradient = torsionAngleGradientRadians(a, b, c, d);

    for(int i = 0; i < gradient.size(); i++){
        gradient[i].scale(chemkit::constants::RadiansToDegrees);
    }

    return gradient;
}

/// Returns the gradient of the torsion angle between the atoms \p a,
/// \p b, \p c, and \p d.
inline QVector<Vector> ForceFieldCalculation::torsionAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return torsionAngleGradientRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline QVector<Vector> ForceFieldCalculation::torsionAngleGradientRadians(const Point &a, const Point &b, const Point &c, const Point &d) const
{
    QVector<Vector> gradient(4);

    Float phi = Point::torsionAngleRadians(a, b, c, d);

    Vector ab = b - a;
    Vector ac = c - a;
    Vector bd = d - b;
    Vector cb = b - c;
    Vector cd = d - c;

    Vector m = ab.cross(cb);
    Vector n = cb.cross(cd);

    Vector p = ((n / (m.length() * n.length())) - ((m / m.lengthSquared()) * cos(phi)));
    Vector q = ((m / (m.length() * n.length())) - ((n / n.lengthSquared()) * cos(phi)));

    gradient[0] = cb.cross(p) * (1.0 / sin(phi));
    gradient[1] = (ac.cross(p) - cd.cross(q)) * (1.0 / sin(phi));
    gradient[2] = (bd.cross(q) - ab.cross(p)) * (1.0 / sin(phi));
    gradient[3] = cb.cross(q) * (1.0 / sin(phi));

    return gradient;
}

/// Returns the wilson angle between the atoms \p a, \p b, \p c, and
/// \p d. The angle is in degrees.
inline Float ForceFieldCalculation::wilsonAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point::wilsonAngle(a->position(), b->position(), c->position(), d->position());
}

/// Returns the wilson angle between the atoms \p a, \p b, \p c, and
/// \p d. The angle is in radians.
inline Float ForceFieldCalculation::wilsonAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return Point::wilsonAngleRadians(a->position(), b->position(), c->position(), d->position());
}

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDCALCULATION_INLINE_H
