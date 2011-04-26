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

#ifndef CHEMKIT_GENERICQUATERNION_H
#define CHEMKIT_GENERICQUATERNION_H

#include "chemkit.h"

#include "genericpoint.h"
#include "genericvector.h"
#include "staticvector.h"

namespace chemkit {

template<typename T>
class GenericQuaternion : public StaticVector<T, 4>
{
    public:
        // construction and destruction
        GenericQuaternion();
        GenericQuaternion(T x, T y, T z, T r);
        GenericQuaternion(const StaticVector<T, 4> &quaternion);

        // properties
        T r() const;
        GenericPoint<T> toPoint3() const;
        GenericVector<T> toVector3() const;

        // math
        GenericQuaternion<T> multiply(const GenericQuaternion<T> &quaternion) const;
        GenericQuaternion<T> conjugate() const;

        // operators
        GenericQuaternion<T> operator*(const GenericQuaternion<T> &quaternion) const;

        // static methods
        static GenericQuaternion<T> rotation(const GenericVector<T> &axis, T angle);
        static GenericQuaternion<T> rotationRadians(const GenericVector<T> &axis, T angle);
        static GenericPoint<T> rotate(const GenericPoint<T> &point, const GenericVector<T> &axis, T angle);
        static GenericPoint<T> rotateRadians(const GenericPoint<T> &point, const GenericVector<T> &axis, T angle);
        static GenericVector<T> rotate(const GenericVector<T> &vector, const GenericVector<T> &axis, T angle);
        static GenericVector<T> rotateRadians(const GenericVector<T> &vector, const GenericVector<T> &axis, T angle);

    private:
        static GenericQuaternion<T> multiply(const GenericQuaternion<T> &a, const GenericQuaternion<T> &b);
};

} // end chemkit namespace

#include "genericquaternion-inline.h"

#endif // CHEMKIT_GENERICQUATERNION_H
