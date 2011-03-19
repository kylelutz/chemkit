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

#ifndef CHEMKIT_COORDINATES_H
#define CHEMKIT_COORDINATES_H

#include "chemkit.h"

#include <QtCore>

#include "matrix.h"
#include "point3.h"
#include "vector3.h"
#include "staticmatrix.h"

namespace chemkit {

class Atom;
class Molecule;
class Conformer;

class CHEMKIT_EXPORT Coordinates
{
    public:
        // construction and destruction
        Coordinates();
        Coordinates(int size);
        Coordinates(const Molecule *molecule);
        Coordinates(const Conformer *conformer);
        Coordinates(const QList<Atom *> &atoms);
        Coordinates(const QList<Point3> &points);
        Coordinates(const Coordinates &coordinates);
        ~Coordinates();

        // properties
        void setSize(int size);
        int size() const;
        bool isEmpty() const;
        Matrix toMatrix() const;

        // coordinates
        void setPosition(int index, const Point3 &position);
        void setPosition(int index, Float x, Float y, Float z);
        Point3 position(int index) const;
        void setValue(int row, int column, Float value);
        Float value(int row, int column) const;
        void append(const Point3 &position);
        void append(Float x, Float y, Float z);
        void insert(int index, const Point3 &position);
        void insert(int index, Float x, Float y, Float z);
        void remove(int index);

        // geometry
        Float distance(int i, int j) const;
        Float angle(int i, int j, int k) const;
        Float angleRadians(int i, int j, int k) const;
        Float torsionAngle(int i, int j, int k, int l) const;
        Float torsionAngleRadians(int i, int j, int k, int l) const;
        Float wilsonAngle(int i, int j, int k, int l) const;
        Float wilsonAngleRadians(int i, int j, int k, int l) const;
        Point3 center() const;
        Point3 weightedCenter(const QVector<Float> &weights) const;
        void moveBy(const Vector3 &vector);
        void moveBy(Float x, Float y, Float z);
        Matrix distanceMatrix() const;

        // math
        Coordinates add(const Coordinates &coordinates) const;
        Coordinates subtract(const Coordinates &coordinates) const;
        StaticMatrix<Float, 3, 3> multiply(const Coordinates *coordinates) const;

        // operators
        Coordinates operator+(const Coordinates &coordinates) const;
        Coordinates operator-(const Coordinates &coordinates) const;
        Coordinates& operator=(const Coordinates &coordinates);

    private:
        Matrix m_matrix;
};

} // end chemkit namespace

#endif // CHEMKIT_COORDINATES_H
