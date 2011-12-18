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

#include "atom.h"
#include "vector3.h"
#include "geometry.h"
#include "molecule.h"

namespace chemkit {

// === CartesianCoordinates ================================================ //
/// \class CartesianCoordinates cartesiancoordinates.h chemkit/cartesiancoordinates.h
/// \ingroup chemkit
/// \brief The CartesianCoordinates class contains cartesian coordinates.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty coordinate matrix.
CartesianCoordinates::CartesianCoordinates()
    : m_matrix(0, 3)
{
}

/// Creates a new, empty coordinate matrix with space for \p size
/// points.
CartesianCoordinates::CartesianCoordinates(int size)
    : m_matrix(size, 3)
{
    m_matrix.setZero();
}

/// Creates a new coordinate matrix with the coordinates from
/// \p molecule.
CartesianCoordinates::CartesianCoordinates(const Molecule *molecule)
    : m_matrix(molecule->size(), 3)
{
    int size = molecule->size();

    for(int i = 0; i < size; i++){
        Point3 position = molecule->atom(i)->position();
        m_matrix(i, 0) = position.x();
        m_matrix(i, 1) = position.y();
        m_matrix(i, 2) = position.z();
    }
}

/// Creates a new coordinate matrix with the coordinates from
/// \p atoms.
CartesianCoordinates::CartesianCoordinates(const std::vector<Atom *> &atoms)
    : m_matrix(atoms.size(), 3)
{
    unsigned int size = atoms.size();

    for(unsigned int i = 0; i < size; i++){
        Point3 position = atoms[i]->position();
        m_matrix(i, 0) = position.x();
        m_matrix(i, 1) = position.y();
        m_matrix(i, 2) = position.z();
    }
}

/// Creates a new coordinate matrix that contains \p points.
CartesianCoordinates::CartesianCoordinates(const std::vector<Point3> &points)
    : m_matrix(points.size(), 3)
{
    for(unsigned int i = 0; i < points.size(); i++){
        Point3 position = points[i];
        m_matrix(i, 0) = position.x();
        m_matrix(i, 1) = position.y();
        m_matrix(i, 2) = position.z();
    }
}

/// Creates a new coordinate matrix that is a copy of \p coordinates.
CartesianCoordinates::CartesianCoordinates(const CartesianCoordinates &coordinates)
    : m_matrix(coordinates.m_matrix)
{
}

/// Destroys the coordinate matrix.
CartesianCoordinates::~CartesianCoordinates()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the size of the matrix to \p size.
void CartesianCoordinates::setSize(int size)
{
    m_matrix.conservativeResize(size, Eigen::NoChange);
}

/// Returns the number of coordinates in the matrix.
int CartesianCoordinates::size() const
{
    return m_matrix.rows();
}

/// Returns \c true if the matrix is empty.
bool CartesianCoordinates::isEmpty() const
{
    return size() == 0;
}

/// Returns a matrix containing the data in the coordinate matrix.
Matrix CartesianCoordinates::toMatrix() const
{
    return m_matrix;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the position at \p index to \p position.
void CartesianCoordinates::setPosition(int index, const Point3 &position)
{
    m_matrix(index, 0) = position.x();
    m_matrix(index, 1) = position.y();
    m_matrix(index, 2) = position.z();
}

/// Sets the position at \p index to (\p x, \p y, \p z).
void CartesianCoordinates::setPosition(int index, Real x, Real y, Real z)
{
    setPosition(index, Point3(x, y, z));
}

/// Returns the coordinates at \p index.
Point3 CartesianCoordinates::position(int index) const
{
    return Point3(m_matrix(index, 0),
                  m_matrix(index, 1),
                  m_matrix(index, 2));
}

/// Sets the value at \p row and \p column to \p value.
void CartesianCoordinates::setValue(int row, int column, Real value)
{
    m_matrix(row, column) = value;
}

/// Returns the value at \p row and \p column;
Real CartesianCoordinates::value(int row, int column) const
{
    return m_matrix(row, column);
}

/// Appends \p position to the coordinates.
void CartesianCoordinates::append(const Point3 &position)
{
    insert(size(), position);
}

/// Appends the point (\p x, \p y, \p z) to the coordinates.
void CartesianCoordinates::append(Real x, Real y, Real z)
{
    append(Point3(x, y, z));
}

/// Inserts \p position at \p index.
void CartesianCoordinates::insert(int index, const Point3 &position)
{
    // resize to make space for the new position
    if(index >= size()){
        setSize(index + 1);
    }
    else{
        setSize(size() + 1);
    }

    // copy old positions
    for(int i = size() - 1; i > index; i--){
        setPosition(i, this->position(i - 1));
    }

    // set the new position
    setPosition(index, position);
}

/// Inserts the point (\p x, \p y, \p z) at \p index.
void CartesianCoordinates::insert(int index, Real x, Real y, Real z)
{
    insert(index, Point3(x, y, z));
}

/// Removes the position at \p index.
void CartesianCoordinates::remove(int index)
{
    for(int i = index + 1; i < size(); i++){
        setPosition(i - 1, position(i));
    }

    setSize(size() - 1);
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the distance between the points at \p i and \p j. The
/// returned distance is in Angstroms.
Real CartesianCoordinates::distance(int i, int j) const
{
    return chemkit::geometry::distance(position(i), position(j));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in degrees.
Real CartesianCoordinates::angle(int i, int j, int k) const
{
    return chemkit::geometry::angle(position(i), position(j), position(k));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in radians.
Real CartesianCoordinates::angleRadians(int i, int j, int k) const
{
    return chemkit::geometry::angleRadians(position(i), position(j), position(k));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in degrees.
Real CartesianCoordinates::torsionAngle(int i, int j, int k, int l) const
{
    return chemkit::geometry::torsionAngle(position(i), position(j), position(k), position(l));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in radians.
Real CartesianCoordinates::torsionAngleRadians(int i, int j, int k, int l) const
{
    return chemkit::geometry::torsionAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in degrees.
Real CartesianCoordinates::wilsonAngle(int i, int j, int k, int l) const
{
    return chemkit::geometry::wilsonAngle(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in radians.
Real CartesianCoordinates::wilsonAngleRadians(int i, int j, int k, int l) const
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

    // sums for each component
    Real sx = 0;
    Real sy = 0;
    Real sz = 0;

    for(int i = 0; i < size(); i++){
        sx += m_matrix(i, 0);
        sy += m_matrix(i, 1);
        sz += m_matrix(i, 2);
    }

    int n = size();

    return Point3(sx/n, sy/n, sz/n);
}

/// Returns the center of the coordinates after weighting each
/// position with \p weights.
Point3 CartesianCoordinates::weightedCenter(const std::vector<Real> &weights) const
{
    assert(static_cast<unsigned long>(size()) == weights.size());

    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Real sx = 0;
    Real sy = 0;
    Real sz = 0;

    // sum of weights
    Real sw = 0;

    for(int i = 0; i < size(); i++){
        Real weight = weights[i];

        sx += weight * m_matrix(i, 0);
        sy += weight * m_matrix(i, 1);
        sz += weight * m_matrix(i, 2);

        sw += weight;
    }

    int n = sw * size();

    return Point3(sx/n, sy/n, sz/n);
}

/// Moves all of the coordinates by \p vector.
void CartesianCoordinates::moveBy(const Vector3 &vector)
{
    for(int i = 0; i < size(); i++){
        m_matrix(i, 0) += vector.x();
        m_matrix(i, 1) += vector.y();
        m_matrix(i, 2) += vector.z();
    }
}

/// Moves all of the coordinates by (\p x, \p y, \p z).
void CartesianCoordinates::moveBy(Real x, Real y, Real z)
{
    moveBy(Vector3(x, y, z));
}

/// Returns a matrix containing the distances between each pair of
/// points in the coordinates.
Matrix CartesianCoordinates::distanceMatrix() const
{
    Matrix matrix(size(), size());

    for(int i = 0; i < size(); i++){
        // set diagonal entries to zero
        matrix(i, i) = 0;

        for(int j = i + 1; j < size(); j++){
            Real d = distance(i, j);

            matrix(i, j) = d;
            matrix(j, i) = d;
        }
    }

    return matrix;
}

// --- Math ---------------------------------------------------------------- //
/// Returns a new coordinate matrix containing the result of adding
/// the coordinates with \p coordinates.
CartesianCoordinates CartesianCoordinates::add(const CartesianCoordinates &coordinates) const
{
    int size = std::min(this->size(), coordinates.size());

    CartesianCoordinates result(size);

    for(int i = 0; i < size; i++){
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
    int size = std::min(this->size(), coordinates.size());

    CartesianCoordinates result(size);

    for(int i = 0; i < size; i++){
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
    return m_matrix.transpose() * coordinates->m_matrix;
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
    m_matrix = coordinates.m_matrix;

    return *this;
}

} // end chemkit namespace
