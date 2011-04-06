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

#include "internalcoordinates.h"

#include <cassert>

#include "vector3.h"
#include "constants.h"
#include "coordinates.h"
#include "staticmatrix.h"

namespace chemkit {

// === InternalCoordinatesPrivate ========================================== //
class InternalCoordinatesPrivate
{
    public:
        int size;
        int *connections;
        Float *coordinates;
};

// === InternalCoordinates ================================================= //
/// \class InternalCoordinates internalcoordinates.h chemkit/internalcoordinates.h
/// \ingroup chemkit
/// \brief The InternalCoordinates class represents a set of internal
///        coordinates.
///
/// \see Coordinates

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty set of internal coordinates.
InternalCoordinates::InternalCoordinates()
    : d(new InternalCoordinatesPrivate)
{
    d->size = 0;
    d->connections = 0;
    d->coordinates = 0;
}

/// Creates a new internal coordinate set with \p size rows.
InternalCoordinates::InternalCoordinates(int size)
    : d(new InternalCoordinatesPrivate)
{
    d->size = size;
    d->connections = new int[3 * size];
    d->coordinates = new Float[3 * size];
}

/// Creates a new internal coordinates object as a copy of
/// \p coordinates.
InternalCoordinates::InternalCoordinates(const InternalCoordinates &coordinates)
    : d(new InternalCoordinatesPrivate)
{
    d->size = coordinates.d->size;
    d->connections = new int[3 * coordinates.d->size];
    d->coordinates = new Float[3 * coordinates.d->size];

    memcpy(d->connections, coordinates.d->connections, 3 * coordinates.d->size * sizeof(int));
    memcpy(d->coordinates, coordinates.d->coordinates, 3 * coordinates.d->size * sizeof(int));
}

/// Destroys the internal coordinates object.
InternalCoordinates::~InternalCoordinates()
{
    delete[] d->connections;
    delete[] d->coordinates;

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of rows of coordinates.
int InternalCoordinates::size() const
{
    return d->size;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in degrees.
void InternalCoordinates::setCoordinates(int row, Float r, Float theta, Float phi)
{
    assert(row < d->size);

    d->coordinates[row * 3 + 0] = r;
    d->coordinates[row * 3 + 1] = theta;
    d->coordinates[row * 3 + 2] = phi;
}

/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in radians.
void InternalCoordinates::setCoordinatesRadians(int row, Float r, Float theta, Float phi)
{
    assert(row < d->size);

    d->coordinates[row * 3 + 0] = r * chemkit::constants::RadiansToDegrees;
    d->coordinates[row * 3 + 1] = theta * chemkit::constants::RadiansToDegrees;
    d->coordinates[row * 3 + 2] = phi * chemkit::constants::RadiansToDegrees;
}

/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in degrees.
std::vector<Float> InternalCoordinates::coordinates(int row) const
{
    assert(row < d->size);

    std::vector<Float> coordinates(3);

    coordinates[0] = d->coordinates[row * 3 + 0];
    coordinates[1] = d->coordinates[row * 3 + 1];
    coordinates[2] = d->coordinates[row * 3 + 2];

    return coordinates;
}

/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in radians.
std::vector<Float> InternalCoordinates::coordinatesRadians(int row) const
{
    assert(row < d->size);

    std::vector<Float> coordinates = this->coordinates(row);

    for(int i = 0; i < 3; i++){
        coordinates[i] *= chemkit::constants::DegreesToRadians;
    }

    return coordinates;
}

/// Sets the connections for the coordinates at \p row to \p a, \p b
/// and \p c.
void InternalCoordinates::setConnections(int row, int a, int b, int c)
{
    assert(row < d->size);

    d->connections[row * 3 + 0] = a;
    d->connections[row * 3 + 1] = b;
    d->connections[row * 3 + 2] = c;
}

/// Returns the connections for the coordinates at \p row.
std::vector<int> InternalCoordinates::connections(int row) const
{
    assert(row < d->size);

    std::vector<int> connections(3);

    connections[0] = d->connections[row * 3 + 0];
    connections[1] = d->connections[row * 3 + 1];
    connections[2] = d->connections[row * 3 + 2];

    return connections;
}

// --- Conversions --------------------------------------------------------- //
/// Converts the internal coordinates into cartesian coordinates.
///
/// The ownership of the returned coordinates object is passed to
/// the caller.
///
/// \internal
/// This method implements the Natural Extension Reference Frame
/// (NeRF) algorithm presented in [Parsons 2005].
/// \endinternal
Coordinates* InternalCoordinates::toCartesianCoordinates() const
{
    Coordinates *cartesianCoordinates = new Coordinates(d->size);

    // set positions for the first three atoms
    if(d->size >= 0){
        cartesianCoordinates->setPosition(0, Point3(0, 0, 0));

        if(d->size >= 1){
            Float r1 = coordinates(1)[0];
            cartesianCoordinates->setPosition(1, Point3(r1, 0, 0));

            if(d->size >= 2){
                Float r2 = coordinates(2)[0];
                Float theta = coordinates(2)[1];

                Float x = r2 * cos((180.0 - theta) * chemkit::constants::DegreesToRadians);
                Float y = r2 * sin((180.0 - theta) * chemkit::constants::DegreesToRadians);

                cartesianCoordinates->setPosition(2, Point3(r1 + x, y, 0));
            }
        }
    }

    // set positions for the rest of the atoms
    for(int i = 3; i < d->size; i++){
        std::vector<Float> coordinates = this->coordinates(i);
        Float r = coordinates[0];
        Float theta = coordinates[1];
        Float phi = coordinates[2];

        Float sinTheta = sin(theta * chemkit::constants::DegreesToRadians);
        Float cosTheta = cos(theta * chemkit::constants::DegreesToRadians);
        Float sinPhi = sin(phi * chemkit::constants::DegreesToRadians);
        Float cosPhi = cos(phi * chemkit::constants::DegreesToRadians);

        Float x = r * cosTheta;
        Float y = r * cosPhi * sinTheta;
        Float z = r * sinPhi * sinTheta;

        std::vector<int> connections = this->connections(i);

        const Point3 &a = cartesianCoordinates->position(connections[2]);
        const Point3 &b = cartesianCoordinates->position(connections[1]);
        const Point3 &c = cartesianCoordinates->position(connections[0]);

        Vector3 ab = (b - a);
        Vector3 bc = (c - b).normalized();
        Vector3 n = ab.cross(bc).normalized();
        Vector3 ncbc = n.cross(bc);

        StaticMatrix<Float, 3, 3> M;
        M << bc.x(), ncbc.x(), n.x(),
             bc.y(), ncbc.y(), n.y(),
             bc.z(), ncbc.z(), n.z();

        Point3 d = M.multiply(Point3(-x, y, z)) + c;
        cartesianCoordinates->setPosition(i, d);
    }

    return cartesianCoordinates;
}

// --- Operators ----------------------------------------------------------- //
InternalCoordinates& InternalCoordinates::operator=(const InternalCoordinates &coordinates)
{
    if(&coordinates == this){
        return *this;
    }

    delete [] d->connections;
    delete [] d->coordinates;

    d->size = coordinates.d->size;
    d->connections = new int[3 * coordinates.d->size];
    d->coordinates = new Float[3 * coordinates.d->size];

    memcpy(d->connections, coordinates.d->connections, 3 * coordinates.d->size * sizeof(int));
    memcpy(d->coordinates, coordinates.d->coordinates, 3 * coordinates.d->size * sizeof(int));

    return *this;
}

} // end chemkit namespace
