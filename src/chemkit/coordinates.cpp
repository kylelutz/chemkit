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

#include "coordinates.h"

#include <cassert>

#include "vector3.h"
#include "molecule.h"

namespace chemkit {

// === Coordinates ========================================================= //
/// \class Coordinates coordinates.h chemkit/coordinates.h
/// \ingroup chemkit
/// \brief The Coordinates class contains cartesian coordinates.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty coordinate matrix.
Coordinates::Coordinates()
    : m_matrix(0, 3)
{
}

/// Creates a new, empty coordinate matrix with space for \p size
/// points.
Coordinates::Coordinates(int size)
    : m_matrix(size, 3)
{
}

/// Creates a new coordinate matrix with the coordinates from
/// \p molecule.
Coordinates::Coordinates(const Molecule *molecule)
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
/// \p conformer.
Coordinates::Coordinates(const Conformer *conformer)
    : m_matrix(conformer->molecule()->size(), 3)
{
    int size = conformer->molecule()->size();

    for(int i = 0; i < size; i++){
        Point3 position = conformer->position(conformer->molecule()->atom(i));
        m_matrix(i, 0) = position.x();
        m_matrix(i, 1) = position.y();
        m_matrix(i, 2) = position.z();
    }
}

/// Creates a new coordinate matrix with the coordinates from
/// \p atoms.
Coordinates::Coordinates(const std::vector<Atom *> &atoms)
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
Coordinates::Coordinates(const std::vector<Point3> &points)
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
Coordinates::Coordinates(const Coordinates &coordinates)
    : m_matrix(coordinates.m_matrix)
{
}

/// Destroys the coordinate matrix.
Coordinates::~Coordinates()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the size of the matrix to \p size.
void Coordinates::setSize(int size)
{
    m_matrix.conservativeResize(size, Eigen::NoChange);
}

/// Returns the number of coordinates in the matrix.
int Coordinates::size() const
{
    return m_matrix.rows();
}

/// Returns \c true if the matrix is empty.
bool Coordinates::isEmpty() const
{
    return size() == 0;
}

/// Returns a matrix containing the data in the coordinate matrix.
Matrix Coordinates::toMatrix() const
{
    return m_matrix;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the position at \p index to \p position.
void Coordinates::setPosition(int index, const Point3 &position)
{
    m_matrix(index, 0) = position.x();
    m_matrix(index, 1) = position.y();
    m_matrix(index, 2) = position.z();
}

/// Sets the position at \p index to (\p x, \p y, \p z).
void Coordinates::setPosition(int index, Float x, Float y, Float z)
{
    setPosition(index, Point3(x, y, z));
}

/// Returns the coordinates at \p index.
Point3 Coordinates::position(int index) const
{
    return Point3(m_matrix(index, 0),
                  m_matrix(index, 1),
                  m_matrix(index, 2));
}

/// Sets the value at \p row and \p column to \p value.
void Coordinates::setValue(int row, int column, Float value)
{
    m_matrix(row, column) = value;
}

/// Returns the value at \p row and \p column;
Float Coordinates::value(int row, int column) const
{
    return m_matrix(row, column);
}

/// Appends \p position to the coordinates.
void Coordinates::append(const Point3 &position)
{
    insert(size(), position);
}

/// Appends the point (\p x, \p y, \p z) to the coordinates.
void Coordinates::append(Float x, Float y, Float z)
{
    append(Point3(x, y, z));
}

/// Inserts \p position at \p index.
void Coordinates::insert(int index, const Point3 &position)
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
void Coordinates::insert(int index, Float x, Float y, Float z)
{
    insert(index, Point3(x, y, z));
}

/// Removes the position at \p index.
void Coordinates::remove(int index)
{
    for(int i = index + 1; i < size(); i++){
        setPosition(i - 1, position(i));
    }

    setSize(size() - 1);
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the distance between the points at \p i and \p j. The
/// returned distance is in Angstroms.
Float Coordinates::distance(int i, int j) const
{
    return Point3::distance(position(i), position(j));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in degrees.
Float Coordinates::angle(int i, int j, int k) const
{
    return Point3::angle(position(i), position(j), position(k));
}

/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in radians.
Float Coordinates::angleRadians(int i, int j, int k) const
{
    return Point3::angleRadians(position(i), position(j), position(k));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in degrees.
Float Coordinates::torsionAngle(int i, int j, int k, int l) const
{
    return Point3::torsionAngle(position(i), position(j), position(k), position(l));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in radians.
Float Coordinates::torsionAngleRadians(int i, int j, int k, int l) const
{
    return Point3::torsionAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in degrees.
Float Coordinates::wilsonAngle(int i, int j, int k, int l) const
{
    return Point3::wilsonAngle(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in radians.
Float Coordinates::wilsonAngleRadians(int i, int j, int k, int l) const
{
    return Point3::wilsonAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the center of the positions in the coordinates. This is
/// also known as the centroid.
Point3 Coordinates::center() const
{
    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Float sx = 0;
    Float sy = 0;
    Float sz = 0;

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
Point3 Coordinates::weightedCenter(const std::vector<Float> &weights) const
{
    assert(static_cast<unsigned long>(size()) == weights.size());

    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Float sx = 0;
    Float sy = 0;
    Float sz = 0;

    // sum of weights
    Float sw = 0;

    for(int i = 0; i < size(); i++){
        Float weight = weights[i];

        sx += weight * m_matrix(i, 0);
        sy += weight * m_matrix(i, 1);
        sz += weight * m_matrix(i, 2);

        sw += weight;
    }

    int n = sw * size();

    return Point3(sx/n, sy/n, sz/n);
}

/// Moves all of the coordinates by \p vector.
void Coordinates::moveBy(const Vector3 &vector)
{
    for(int i = 0; i < size(); i++){
        m_matrix(i, 0) += vector.x();
        m_matrix(i, 1) += vector.y();
        m_matrix(i, 2) += vector.z();
    }
}

/// Moves all of the coordinates by (\p x, \p y, \p z).
void Coordinates::moveBy(Float x, Float y, Float z)
{
    moveBy(Vector3(x, y, z));
}

/// Returns a matrix containing the distances between each pair of
/// points in the coordinates.
Matrix Coordinates::distanceMatrix() const
{
    Matrix matrix(size(), size());

    for(int i = 0; i < size(); i++){
        // set diagonal entries to zero
        matrix(i, i) = 0;

        for(int j = i + 1; j < size(); j++){
            Float d = distance(i, j);

            matrix(i, j) = d;
            matrix(j, i) = d;
        }
    }

    return matrix;
}

// --- Math ---------------------------------------------------------------- //
/// Returns a new coordinate matrix containing the result of adding
/// the coordinates with \p coordinates.
Coordinates Coordinates::add(const Coordinates &coordinates) const
{
    int size = std::min(this->size(), coordinates.size());

    Coordinates result(size);

    for(int i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a + b);
    }

    return result;
}

/// Returns a new coordinate matrix containing the result of
/// subtracting the coordinates with \p coordinates.
Coordinates Coordinates::subtract(const Coordinates &coordinates) const
{
    int size = std::min(this->size(), coordinates.size());

    Coordinates result(size);

    for(int i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a - b);
    }

    return result;
}

/// Returns the 3x3 matrix product of the transpose of the matrix
/// and \p coordinates.
StaticMatrix<Float, 3, 3> Coordinates::multiply(const Coordinates *coordinates) const
{
    Eigen::Matrix<Float, 3, 3> p = m_matrix.transpose() * coordinates->m_matrix;

    StaticMatrix<Float, 3, 3> product;

    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            product(i, j) = p(i, j);
        }
    }

    return product;
}

// --- Operators ----------------------------------------------------------- //
Coordinates Coordinates::operator+(const Coordinates &coordinates) const
{
    return add(coordinates);
}

Coordinates Coordinates::operator-(const Coordinates &coordinates) const
{
    return subtract(coordinates);
}

Coordinates& Coordinates::operator=(const Coordinates &coordinates)
{
    m_matrix = coordinates.m_matrix;

    return *this;
}

} // end chemkit namespace
