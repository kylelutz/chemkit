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
inline boost::array<Vector3, 2> ForceFieldCalculation::distanceGradient(const ForceFieldAtom *a, const ForceFieldAtom *b) const
{
    return chemkit::geometry::distanceGradient(a->position(), b->position());
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
inline boost::array<Vector3, 3> ForceFieldCalculation::bondAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return chemkit::geometry::angleGradient(a->position(), b->position(), c->position());
}

/// Returns the gradient of the bond angle between atoms \p a, \p b
/// and \p c.
inline boost::array<Vector3, 3> ForceFieldCalculation::bondAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c) const
{
    return chemkit::geometry::angleGradientRadians(a->position(), b->position(), c->position());
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
inline boost::array<Vector3, 4> ForceFieldCalculation::torsionAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::torsionAngleGradient(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the torsion angle between the atoms \p a,
/// \p b, \p c, and \p d.
inline boost::array<Vector3, 4> ForceFieldCalculation::torsionAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::torsionAngleGradientRadians(a->position(), b->position(), c->position(), d->position());
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
inline boost::array<Vector3, 4> ForceFieldCalculation::wilsonAngleGradient(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::wilsonAngleGradient(a->position(), b->position(), c->position(), d->position());
}

/// Returns the gradient of the wilson angle between the atoms
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> ForceFieldCalculation::wilsonAngleGradientRadians(const ForceFieldAtom *a, const ForceFieldAtom *b, const ForceFieldAtom *c, const ForceFieldAtom *d) const
{
    return chemkit::geometry::wilsonAngleGradientRadians(a->position(), b->position(), c->position(), d->position());
}

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDCALCULATION_INLINE_H
