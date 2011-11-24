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

#ifndef CHEMKIT_FORCEFIELDCALCULATION_INLINE_H
#define CHEMKIT_FORCEFIELDCALCULATION_INLINE_H

#include "forcefieldcalculation.h"

#include <chemkit/point3.h>
#include <chemkit/geometry.h>

#include "forcefieldatom.h"

namespace chemkit {

/// Returns the distance between atoms \p a and \p b. Distance is in
/// Angstroms.
inline Real ForceFieldCalculation::distance(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return chemkit::geometry::distance(a->position(), b->position());
}

/// Returns the gradient of the distance between atoms \p a and \p b.
inline std::vector<Vector3> ForceFieldCalculation::distanceGradient(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return distanceGradient(a->position(), b->position());
}

/// Returns the gradient of the distance between points \p a and \p b.
inline std::vector<Vector3> ForceFieldCalculation::distanceGradient(const Point3 &a, const Point3 &b) const
{
    std::vector<Vector3> gradient(2);

    Real distance = chemkit::geometry::distance(a, b);

    gradient[0] = (a - b) / distance;
    gradient[1] = -gradient[0];

    return gradient;
}

/// Returns the bond angle between atoms \p a, \p b and \p c. The
/// angle is in degrees.
inline Real ForceFieldCalculation::bondAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return chemkit::geometry::angle(a->position(), b->position(), c->position());
}

/// Returns the bond angle between atoms \p a, \p b and \p c. The
/// angle is in radians.
inline Real ForceFieldCalculation::bondAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return chemkit::geometry::angleRadians(a->position(), b->position(), c->position());
}

/// Returns the gradient of the bond angle between atoms \p a, \p b
/// and \p c.
inline std::vector<Vector3> ForceFieldCalculation::bondAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    std::vector<Vector3> gradient = bondAngleGradientRadians(a, b, c);

    for(unsigned int i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the bond angle between atoms \p a, \p b
/// and \p c.
inline std::vector<Vector3> ForceFieldCalculation::bondAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return bondAngleGradientRadians(a->position(), b->position(), c->position());
}

/// Returns the gradient of the bond angle between points \p a, \p b
/// and \p c.
inline std::vector<Vector3> ForceFieldCalculation::bondAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c) const
{
    std::vector<Vector3> gradient(3);

    Real theta = chemkit::geometry::angleRadians(a, b, c);

    Real rab = chemkit::geometry::distance(a, b);
    Real rbc = chemkit::geometry::distance(b, c);

    gradient[0] = ((((c - b) * rab) - (a - b) * ((b - a).dot(b - c) / rab)) / (pow(rab, 2) * rbc)) / -sin(theta);
    gradient[1] = ((((b - c) + (b - a)) * (rab * rbc) - (((b - a) * (rbc/rab) + (b - c) * (rab/rbc)) * (b - a).dot(b - c))) / pow(rab * rbc, 2)) / -sin(theta);
    gradient[2] = -gradient[0] - gradient[1];

    return gradient;
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between the atoms \p a, \p b, \p c, and \p d. The angle is in
/// degrees.
inline Real ForceFieldCalculation::torsionAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::torsionAngle(a->position(), b->position(), c->position(), d->position());
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between the atoms \p a, \p b, \p c, and \p d. The angle is in
/// radians.
inline Real ForceFieldCalculation::torsionAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::torsionAngleRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the torsion angle between the atoms \p a,
/// \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::torsionAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    std::vector<Vector3> gradient = torsionAngleGradientRadians(a, b, c, d);

    for(unsigned int i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the torsion angle between the atoms \p a,
/// \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::torsionAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return torsionAngleGradientRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::torsionAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d) const
{
    std::vector<Vector3> gradient(4);

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

/// Returns the wilson angle between the atoms \p a, \p b, \p c, and
/// \p d. The angle is in degrees.
inline Real ForceFieldCalculation::wilsonAngle(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::wilsonAngle(a->position(), b->position(), c->position(), d->position());
}

/// Returns the wilson angle between the atoms \p a, \p b, \p c, and
/// \p d. The angle is in radians.
inline Real ForceFieldCalculation::wilsonAngleRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::wilsonAngleRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the wilson angle between the atoms
/// \p a, \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::wilsonAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    std::vector<Vector3> gradient = wilsonAngleGradientRadians(a, b, c, d);

    for(unsigned int i = 0; i < gradient.size(); i++){
        gradient[i] *= chemkit::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the wilson angle between the atoms
/// \p a, \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::wilsonAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return wilsonAngleGradientRadians(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the wilson angle between the points
/// \p a, \p b, \p c, and \p d.
inline std::vector<Vector3> ForceFieldCalculation::wilsonAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d) const
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

    std::vector<Vector3> gradient(4);

    gradient[0] = ((bd.cross(bc) / (cos(w) * sin(theta)) - (ba - bc * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rba;
    gradient[2] = ((ba.cross(bd) / (cos(w) * sin(theta)) - (bc - ba * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rbc;
    gradient[3] = (bc.cross(ba) / (cos(w) * sin(theta)) - bd * tan(w)) / rbd;
    gradient[1] = -(gradient[0] + gradient[2] + gradient[3]);

    return gradient;
}

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDCALCULATION_INLINE_H
