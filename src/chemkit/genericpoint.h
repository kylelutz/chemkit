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

#ifndef CHEMKIT_GENERICPOINT_H
#define CHEMKIT_GENERICPOINT_H

#include "staticvector.h"

namespace chemkit {

template<typename T>
class GenericPoint : public StaticVector<T, 3>
{
    public:
        // construction and destruction
        GenericPoint();
        GenericPoint(T x, T y, T z);
        GenericPoint(const StaticVector<float, 3> &vector);
        GenericPoint(const StaticVector<double, 3> &vector);

        // properties
        void moveBy(T dx, T dy, T dz);
        void moveBy(const StaticVector<T, 3> &vector);
        void moveBy(T distance, const StaticVector<T, 3> &vector);
        GenericPoint<T> movedBy(T dx, T dy, T dz) const;
        GenericPoint<T> movedBy(const StaticVector<T, 3> &vector) const;
        GenericPoint<T> movedBy(T distance, const StaticVector<T, 3> &direction) const;

        // math
        T distance(const GenericPoint<T> &point) const;
        GenericPoint<T> midpoint(const GenericPoint<T> &point) const;

        // static methods
        static T angle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c);
        static T angleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c);
        static T torsionAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T torsionAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T wilsonAngle(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static T wilsonAngleRadians(const GenericPoint<T> &a, const GenericPoint<T> &b, const GenericPoint<T> &c, const GenericPoint<T> &d);
        static GenericPoint<T> midpoint(const GenericPoint<T> &a, const GenericPoint<T> &b);
};

} // end chemkit namespace

#include "genericpoint-inline.h"

#endif // CHEMKIT_GENERICPOINT_H
