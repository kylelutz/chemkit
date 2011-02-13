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

#ifndef CHEMKIT_GENERICVECTOR_H
#define CHEMKIT_GENERICVECTOR_H

#include "staticvector.h"

namespace chemkit {

template<typename T>
class GenericVector : public StaticVector<T, 3>
{
    public:
        // construction and destruction
        GenericVector();
        GenericVector(T x, T y, T z);
        GenericVector(const StaticVector<float, 3> &vector);
        GenericVector(const StaticVector<double, 3> &vector);

        // properties
        T x() const;
        T y() const;
        T z() const;

        // static methods
        static GenericVector<T> X();
        static GenericVector<T> Y();
        static GenericVector<T> Z();
        static GenericVector<T> planeNormal(const StaticVector<T, 3> &a, const StaticVector<T, 3> &b, const StaticVector<T, 3> &c);
        static GenericVector<T> randomUnitVector();
};

} // end chemkit namespace

#include "genericvector-inline.h"

#endif // CHEMKIT_GENERICVECTOR_H
