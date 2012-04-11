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

#include "cartesiancoordinates.h"

#include <cassert>

#include "vector3.h"
#include "foreach.h"
#include "geometry.h"

namespace chemkit {

// === CartesianCoordinates ================================================ //
/// \class CartesianCoordinates cartesiancoordinates.h chemkit/cartesiancoordinates.h
/// \ingroup chemkit
/// \brief The CartesianCoordinates class contains cartesian coordinates.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty coordinate matrix.
CartesianCoordinates::CartesianCoordinates()
{
}

/// Creates a new, empty coordinate matrix with space for \p size
/// points.
CartesianCoordinates::CartesianCoordinates(size_t size)
    : m_coordinates(size)
{
    std::fill(m_coordinates.begin(), m_coordinates.end(), Point3(0, 0, 0));
}

/// Creates a new coordinate matrix that is a copy of \p coordinates.
CartesianCoordinates::CartesianCoordinates(const CartesianCoordinates &coordinates)
    : m_coordinates(coordinates.m_coordinates)
{
}

/// Destroys the coordinate matrix.
CartesianCoordinates::~CartesianCoordinates()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the size of the matrix to \p size.
void CartesianCoordinates::resize(size_t size)
{
    m_coordinates.resize(size);
}

/// Returns the number of coordinates in the matrix.
size_t CartesianCoordinates::size() const
{
    return m_coordinates.size();
}

/// Returns \c true if the matrix is empty.
bool CartesianCoordinates::isEmpty() const
{
    return m_coordinates.empty();
}

/// Returns a matrix containing the data in the coordinate matrix.
Matrix CartesianCoordinates::toMatrix() const
{
    Matrix matrix(size(), 3);

    for(size_t i = 0; i < m_coordinates.size(); i++){
        const Point3 &point = m_coordinates[i];

        matrix(i, 0) = point.x();
        matrix(i, 1) = point.y();
        matrix(i, 2) = point.z();
    }

    return matrix;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the position at \p index to \p position.
void CartesianCoordinates::setPosition(size_t index, const Point3 &position)
{
    assert(index < size());

    m_coordinates[index] = position;
}

/// Sets the position at \p index to (\p x, \p y, \p z).
void CartesianCoordinates::setPosition(size_t index, Real x, Real y, Real z)
{
    setPosition(index, Point3(x, y, z));
}

/// Returns the coordinates at \p index.
Point3 CartesianCoordinates::position(size_t index) const
{
    assert(index < size());

    return m_coordinates[index];
}

/// Sets the value at \p row and \p column to \p value.
void CartesianCoordinates::setValue(size_t row, size_t column, Real value)
{
    assert(row < size());
    assert(column < 3);

    m_coordinates[row][column] = value;
}

/// Returns the value at \p row and \p column;
Real CartesianCoordinates::value(size_t row, size_t column) const
{
    assert(row < size());
    assert(column < 3);

    return m_coordinates[row][column];
}

/// Appends \p position to the coordinates.
void CartesianCoordinates::append(const Point3 &position)
{
    m_coordinates.push_back(position);
}

/// Appends the point (\p x, \p y, \p z) to the coordinates.
void CartesianCoordinates::append(Real x, Real y, Real z)
{
    append(Point3(x, y, z));
}

/// Inserts \p position at \p index.
void CartesianCoordinates::insert(size_t index, const Point3 &position)
{
    m_coordinates.insert(m_coordinates.begin() + index, position);
}

/// Inserts the point (\p x, \p y, \p z) at \p index.
void CartesianCoordinates::insert(size_t index, Real x, Real y, Real z)
{
    insert(index, Point3(x, y, z));
}

/// Removes the position at \p index.
void CartesianCoordinates::remove(size_t index)
{
    m_coordinates.erase(m_coordinates.begin() + index);
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the distance between the points at \p i and \p j. The
/// returned distance is in Angstroms.
Real CartesianCoordinates::distance(size_t i, size_t j) const
{
    return chemkit::geometry::distance(position(i), position(j));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in degrees.
Real CartesianCoordinates::angle(size_t i, size_t j, size_t k) const
{
    return chemkit::geometry::angle(position(i), position(j), position(k));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in radians.
Real CartesianCoordinates::angleRadians(size_t i, size_t j, size_t k) const
{
    return chemkit::geometry::angleRadians(position(i), position(j), position(k));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in degrees.
Real CartesianCoordinates::torsionAngle(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::torsionAngle(position(i), position(j), position(k), position(l));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in radians.
Real CartesianCoordinates::torsionAngleRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::torsionAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in degrees.
Real CartesianCoordinates::wilsonAngle(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::wilsonAngle(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in radians.
Real CartesianCoordinates::wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::wilsonAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the center of the positions in the coordinates. This is
/// also known as the centroid.
Point3 CartesianCoordinates::center() const
{
    if(isEmpty()){
        return Point3(0, 0, 0);
    }

    Point3 sum(0, 0, 0);

    for(size_t i = 0; i < size(); i++){
        sum += m_coordinates[i];
    }

    return (1.0 / size()) * sum;
}

/// Returns the center of the coordinates after weighting each
/// position with \p weights.
Point3 CartesianCoordinates::weightedCenter(const std::vector<Real> &weights) const
{
    assert(weights.size() == size());

    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Point3 sum(0, 0, 0);

    // sum of weights
    Real sw = 0;

    for(size_t i = 0; i < size(); i++){
        sum += weights[i] * m_coordinates[i];
        sw += weights[i];
    }

    return (1.0 / (sw * size())) * sum;
}

/// Moves all of the coordinates by \p vector.
void CartesianCoordinates::moveBy(const Vector3 &vector)
{
    foreach(Point3 &point, m_coordinates){
        point += vector;
    }
}

/// Moves all of the coordinates by (\p x, \p y, \p z).
void CartesianCoordinates::moveBy(Real x, Real y, Real z)
{
    moveBy(Vector3(x, y, z));
}

/// Rotates the coordinates by \p angle degrees around \p axis.
void CartesianCoordinates::rotate(const Vector3 &axis, Real angle)
{
    // convert angle to radians
    angle *= chemkit::constants::DegreesToRadians;

    // build rotation transform
    Eigen::Matrix<Real, 3, 1> axisVector(axis.x(), axis.y(), axis.z());
    Eigen::Transform<Real, 3, 3> transform(Eigen::AngleAxis<Real>(angle, axisVector));

    // rotate each point
    for(size_t i = 0; i < m_coordinates.size(); i++){
        setPosition(i, transform * position(i));
    }
}

/// Returns a matrix containing the distances between each pair of
/// points in the coordinates.
Matrix CartesianCoordinates::distanceMatrix() const
{
    Matrix matrix(size(), size());

    for(size_t i = 0; i < size(); i++){
        // set diagonal entries to zero
        matrix(i, i) = 0;

        for(size_t j = i + 1; j < size(); j++){
            Real d = distance(i, j);

            matrix(i, j) = d;
            matrix(j, i) = d;
        }
    }

    return matrix;
}

// --- Derivatives --------------------------------------------------------- //
/// Returns the gradient of the distance between the points at \p i
/// and \p j.
boost::array<Vector3, 2> CartesianCoordinates::distanceGradient(size_t i, size_t j) const
{
    return chemkit::geometry::distanceGradient(position(i), position(j));
}

/// Returns the gradient of the angle between the points at \p i,
/// \p j and \p k.
boost::array<Vector3, 3> CartesianCoordinates::angleGradient(size_t i, size_t j, size_t k) const
{
    return chemkit::geometry::angleGradient(position(i),
                                            position(j),
                                            position(k));
}

/// Returns the gradient of the angle between the points at \p i,
/// \p j and \p k.
boost::array<Vector3, 3> CartesianCoordinates::angleGradientRadians(size_t i, size_t j, size_t k) const
{
    return chemkit::geometry::angleGradientRadians(position(i),
                                                   position(j),
                                                   position(k));
}

/// Returns the gradient of the torsion angle between the points
/// at \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CartesianCoordinates::torsionAngleGradient(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::torsionAngleGradient(position(i),
                                                   position(j),
                                                   position(k),
                                                   position(l));
}

/// Returns the gradient of the torsion angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CartesianCoordinates::torsionAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::torsionAngleGradientRadians(position(i),
                                                          position(j),
                                                          position(k),
                                                          position(l));
}

/// Returns the gradient of the wilson angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CartesianCoordinates::wilsonAngleGradient(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::wilsonAngleGradient(position(i),
                                                  position(j),
                                                  position(k),
                                                  position(l));
}

/// Returns the gradient of the wilson angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CartesianCoordinates::wilsonAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return chemkit::geometry::wilsonAngleGradientRadians(position(i),
                                                         position(j),
                                                         position(k),
                                                         position(l));
}

// --- Math ---------------------------------------------------------------- //
/// Returns a new coordinate matrix containing the result of adding
/// the coordinates with \p coordinates.
CartesianCoordinates CartesianCoordinates::add(const CartesianCoordinates &coordinates) const
{
    size_t size = std::min(this->size(), coordinates.size());

    CartesianCoordinates result(size);

    for(size_t i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a + b);
    }

    return result;
}

/// Returns a new coordinate matrix containing the result of
/// subtracting the coordinates with \p coordinates.
CartesianCoordinates CartesianCoordinates::subtract(const CartesianCoordinates &coordinates) const
{
    size_t size = std::min(this->size(), coordinates.size());

    CartesianCoordinates result(size);

    for(size_t i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a - b);
    }

    return result;
}

/// Returns the 3x3 matrix product of the transpose of the matrix
/// and \p coordinates.
Eigen::Matrix<Real, 3, 3> CartesianCoordinates::multiply(const CartesianCoordinates *coordinates) const
{
    assert(coordinates->size() == this->size());

    return toMatrix().transpose() * coordinates->toMatrix();
}

// --- Operators ----------------------------------------------------------- //
CartesianCoordinates CartesianCoordinates::operator+(const CartesianCoordinates &coordinates) const
{
    return add(coordinates);
}

CartesianCoordinates CartesianCoordinates::operator-(const CartesianCoordinates &coordinates) const
{
    return subtract(coordinates);
}

CartesianCoordinates& CartesianCoordinates::operator=(const CartesianCoordinates &coordinates)
{
    if(this != &coordinates){
        m_coordinates = coordinates.m_coordinates;
    }

    return *this;
}

/// Returns the position at \p index.
Point3& CartesianCoordinates::operator[](size_t index)
{
    return m_coordinates[index];
}

/// \overload
const Point3& CartesianCoordinates::operator[](size_t index) const
{
    return m_coordinates[index];
}

} // end chemkit namespace
