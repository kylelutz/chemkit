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

namespace chemkit {

namespace geometry {

// --- Transforms ---------------------------------------------------------- //
template<typename T>
inline StaticVector<T, 3> rotate(const StaticVector<T, 3> &vector, const StaticVector<T, 3> &axis, T angle)
{
    return rotateRadians<T>(vector, axis, angle * chemkit::constants::DegreesToRadians);
}

template<typename T>
inline StaticVector<T, 3> rotateRadians(const StaticVector<T, 3> &vector, const StaticVector<T, 3> &axis, T angle)
{
    Eigen::Matrix<T, 3, 1> axisVector(axis.x(), axis.y(), axis.z());
    Eigen::Transform<T, 3, 3> transform(Eigen::AngleAxis<T>(angle, axisVector));

    Eigen::Matrix<T, 3, 1> v(vector.x(), vector.y(), vector.z());

    v = transform * v;

    StaticVector<T, 3> result;
    result[0] = v[0];
    result[1] = v[1];
    result[2] = v[2];

    return result;
}

} // end geometry namespace

} // end chemkit namespace

#endif // CHEMKIT_GEOMETRY_INLINE_H
