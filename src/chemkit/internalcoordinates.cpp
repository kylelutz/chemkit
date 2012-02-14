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

#include "internalcoordinates.h"

#include <cassert>

#include "vector3.h"
#include "constants.h"
#include "cartesiancoordinates.h"

namespace chemkit {

// === InternalCoordinatesPrivate ========================================== //
class InternalCoordinatesPrivate
{
public:
    size_t size;
    size_t *connections;
    Real *coordinates;
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
InternalCoordinates::InternalCoordinates(size_t size)
    : d(new InternalCoordinatesPrivate)
{
    d->size = size;
    d->connections = new size_t[3 * size];
    d->coordinates = new Real[3 * size];
}

/// Creates a new internal coordinates object as a copy of
/// \p coordinates.
InternalCoordinates::InternalCoordinates(const InternalCoordinates &coordinates)
    : d(new InternalCoordinatesPrivate)
{
    d->size = coordinates.d->size;
    d->connections = new size_t[3 * coordinates.d->size];
    d->coordinates = new Real[3 * coordinates.d->size];

    memcpy(d->connections, coordinates.d->connections, 3 * coordinates.d->size * sizeof(size_t));
    memcpy(d->coordinates, coordinates.d->coordinates, 3 * coordinates.d->size * sizeof(size_t));
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
size_t InternalCoordinates::size() const
{
    return d->size;
}

/// Returns \c true if the internal coordinates object contains
/// no coordinates (i.e. size() == \c 0).
bool InternalCoordinates::isEmpty() const
{
    return size() == 0;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in degrees.
void InternalCoordinates::setCoordinates(size_t row, Real r, Real theta, Real phi)
{
    assert(row < d->size);

    d->coordinates[row * 3 + 0] = r;
    d->coordinates[row * 3 + 1] = theta;
    d->coordinates[row * 3 + 2] = phi;
}

/// Sets the distance, angle, and torsion at \p row to \p r,
/// \p theta and \p phi respectively. The angles are in radians.
void InternalCoordinates::setCoordinatesRadians(size_t row, Real r, Real theta, Real phi)
{
    assert(row < d->size);

    d->coordinates[row * 3 + 0] = r * chemkit::constants::RadiansToDegrees;
    d->coordinates[row * 3 + 1] = theta * chemkit::constants::RadiansToDegrees;
    d->coordinates[row * 3 + 2] = phi * chemkit::constants::RadiansToDegrees;
}

/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in degrees.
std::vector<Real> InternalCoordinates::coordinates(size_t row) const
{
    assert(row < d->size);

    std::vector<Real> coordinates(3);

    coordinates[0] = d->coordinates[row * 3 + 0];
    coordinates[1] = d->coordinates[row * 3 + 1];
    coordinates[2] = d->coordinates[row * 3 + 2];

    return coordinates;
}

/// Returns the distance, angle, and torsion coordinates at \p row.
/// The returned angles are in radians.
std::vector<Real> InternalCoordinates::coordinatesRadians(size_t row) const
{
    assert(row < d->size);

    std::vector<Real> coordinates = this->coordinates(row);

    for(size_t i = 0; i < 3; i++){
        coordinates[i] *= chemkit::constants::DegreesToRadians;
    }

    return coordinates;
}

/// Sets the connections for the coordinates at \p row to \p a, \p b
/// and \p c.
void InternalCoordinates::setConnections(size_t row, size_t a, size_t b, size_t c)
{
    assert(row < d->size);

    d->connections[row * 3 + 0] = a;
    d->connections[row * 3 + 1] = b;
    d->connections[row * 3 + 2] = c;
}

/// Returns the connections for the coordinates at \p row.
std::vector<size_t> InternalCoordinates::connections(size_t row) const
{
    assert(row < d->size);

    std::vector<size_t> connections(3);

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
/// This method implements the
/// \blueobeliskalgorithm{zmatrixCoordinatesIntoCartesianCoordinates}.
///
/// \internal
/// This method implements the Natural Extension Reference Frame
/// (NeRF) algorithm presented in [Parsons 2005].
/// \endinternal
CartesianCoordinates* InternalCoordinates::toCartesianCoordinates() const
{
    CartesianCoordinates *cartesianCoordinates = new CartesianCoordinates(d->size);

    // set positions for the first three atoms
    if(d->size > 0){
        cartesianCoordinates->setPosition(0, Point3(0, 0, 0));

        if(d->size > 1){
            Real r1 = coordinates(1)[0];
            cartesianCoordinates->setPosition(1, Point3(r1, 0, 0));

            if(d->size > 2){
                Real r2 = coordinates(2)[0];
                Real theta = coordinates(2)[1];

                Real x = r2 * cos((180.0 - theta) * chemkit::constants::DegreesToRadians);
                Real y = r2 * sin((180.0 - theta) * chemkit::constants::DegreesToRadians);

                cartesianCoordinates->setPosition(2, Point3(r1 + x, y, 0));
            }
        }
    }

    // set positions for the rest of the atoms
    for(size_t i = 3; i < d->size; i++){
        std::vector<Real> coordinates = this->coordinates(i);
        Real r = coordinates[0];
        Real theta = coordinates[1];
        Real phi = coordinates[2];

        Real sinTheta = sin(theta * chemkit::constants::DegreesToRadians);
        Real cosTheta = cos(theta * chemkit::constants::DegreesToRadians);
        Real sinPhi = sin(phi * chemkit::constants::DegreesToRadians);
        Real cosPhi = cos(phi * chemkit::constants::DegreesToRadians);

        Real x = r * cosTheta;
        Real y = r * cosPhi * sinTheta;
        Real z = r * sinPhi * sinTheta;

        std::vector<size_t> connections = this->connections(i);

        const Point3 &a = cartesianCoordinates->position(connections[2]);
        const Point3 &b = cartesianCoordinates->position(connections[1]);
        const Point3 &c = cartesianCoordinates->position(connections[0]);

        Vector3 ab = (b - a);
        Vector3 bc = (c - b).normalized();
        Vector3 n = ab.cross(bc).normalized();
        Vector3 ncbc = n.cross(bc);

        Eigen::Matrix<Real, 3, 3> M;
        M << bc.x(), ncbc.x(), n.x(),
             bc.y(), ncbc.y(), n.y(),
             bc.z(), ncbc.z(), n.z();

        Point3 d = (M * Point3(-x, y, z)) + c;
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
    d->connections = new size_t[3 * coordinates.d->size];
    d->coordinates = new Real[3 * coordinates.d->size];

    memcpy(d->connections, coordinates.d->connections, 3 * coordinates.d->size * sizeof(size_t));
    memcpy(d->coordinates, coordinates.d->coordinates, 3 * coordinates.d->size * sizeof(size_t));

    return *this;
}

} // end chemkit namespace
