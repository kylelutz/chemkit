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

#ifndef CHEMKIT_SCALARFIELD_H
#define CHEMKIT_SCALARFIELD_H

#include "chemkit.h"

#include <QtCore>

#include "point3.h"
#include "vector.h"

namespace chemkit {

class ScalarFieldPrivate;

class CHEMKIT_EXPORT ScalarField
{
    public:
        // construction and destruction
        ScalarField();
        ScalarField(const QVector<int> &dimensions, const QVector<Float> &cellLengths, const QVector<Float> &data);
        ~ScalarField();

        // properties
        int width() const;
        int height() const;
        int depth() const;
        int size() const;
        QVector<int> dimensions() const;
        Float cellWidth() const;
        Float cellHeight() const;
        Float cellDepth() const;
        QVector<Float> cellDimensions() const;
        void setOrigin(const Point3 &origin);
        Point3 origin() const;
        QVector<Float> data() const;

        // values
        void setValue(int i, int j, int k, Float value);
        Float value(int i, int j, int k) const;
        Float value(const Point3 &position) const;
        Point3 position(int i, int j, int k) const;
        Vector gradient(int i, int j, int k) const;
        Vector gradient(const Point3 &position) const;

    private:
        ScalarFieldPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_SCALARFIELD_H
