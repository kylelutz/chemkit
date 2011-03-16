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

#ifndef CHEMKIT_INTERNALCOORDINATES_H
#define CHEMKIT_INTERNALCOORDINATES_H

#include "chemkit.h"

#include "point.h"

namespace chemkit {

class Coordinates;
class InternalCoordinatesPrivate;

class CHEMKIT_EXPORT InternalCoordinates
{
    public:
        // construction and destruction
        InternalCoordinates();
        InternalCoordinates(int size);
        InternalCoordinates(const InternalCoordinates &coordinates);
        ~InternalCoordinates();

        // properties
        int size() const;

        // coordinates
        void setCoordinates(int row, Float r, Float theta = 0, Float phi = 0);
        void setCoordinatesRadians(int row, Float r, Float theta = 0, Float phi = 0);
        QVector<Float> coordinates(int row) const;
        QVector<Float> coordinatesRadians(int row) const;
        void setConnections(int row, int a, int b = 0, int c = 0);
        QVector<int> connections(int row) const;

        // conversions
        Coordinates* toCartesianCoordinates() const;

        // operators
        InternalCoordinates& operator=(const InternalCoordinates &coordinates);

    private:
        InternalCoordinatesPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_INTERNALCOORDINATES_H
