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

#ifndef CHEMKIT_COORDINATEMATRIX_H
#define CHEMKIT_COORDINATEMATRIX_H

#include "chemkit.h"

#include <QtCore>

#include "point.h"
#include "matrix.h"
#include "staticmatrix.h"

namespace chemkit {

class Atom;
class Vector;
class Molecule;
class Conformer;

class CHEMKIT_EXPORT CoordinateMatrix
{
    public:
        // construction and destruction
        CoordinateMatrix();
        CoordinateMatrix(int size);
        CoordinateMatrix(const Molecule *molecule);
        CoordinateMatrix(const Conformer *conformer);
        CoordinateMatrix(const QList<Atom *> &atoms);
        CoordinateMatrix(const QList<Point> &points);
        CoordinateMatrix(const CoordinateMatrix &matrix);
        ~CoordinateMatrix();

        // properties
        void setSize(int size);
        int size() const;
        bool isEmpty() const;
        Matrix toMatrix() const;

        // coordinates
        void setPosition(int index, const Point &position);
        void setPosition(int index, Float x, Float y, Float z);
        Point position(int index) const;
        void setValue(int row, int column, Float value);
        Float value(int row, int column) const;
        void append(const Point &position);
        void append(Float x, Float y, Float z);
        void insert(int index, const Point &position);
        void insert(int index, Float x, Float y, Float z);
        void remove(int index);

        // geometry
        Float distance(int i, int j) const;
        Float bondAngle(int i, int j, int k) const;
        Float torsionAngle(int i, int j, int k, int l) const;
        Float wilsonAngle(int i, int j, int k, int l) const;
        Point center() const;
        Point weightedCenter(const QVector<Float> &weights) const;
        void moveBy(const Vector &vector);
        void moveBy(Float x, Float y, Float z);
        Matrix distanceMatrix() const;

        // math
        CoordinateMatrix add(const CoordinateMatrix &matrix) const;
        CoordinateMatrix subtract(const CoordinateMatrix &matrix) const;
        StaticMatrix<Float, 3, 3> multiply(const CoordinateMatrix *matrix) const;

        // operators
        CoordinateMatrix operator+(const CoordinateMatrix &matrix) const;
        CoordinateMatrix operator-(const CoordinateMatrix &matrix) const;
        CoordinateMatrix& operator=(const CoordinateMatrix &matrix);

    private:
        Matrix m_matrix;
};

} // end chemkit namespace

#endif // CHEMKIT_COORDINATEMATRIX_H
