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

#ifndef CHEMKIT_CARTESIANCOORDINATES_H
#define CHEMKIT_CARTESIANCOORDINATES_H

#include "chemkit.h"

#include <vector>

#include <boost/array.hpp>

#include "matrix.h"
#include "point3.h"
#include "vector3.h"

namespace chemkit {

class CHEMKIT_EXPORT CartesianCoordinates
{
public:
    // construction and destruction
    CartesianCoordinates();
    CartesianCoordinates(size_t size);
    CartesianCoordinates(const CartesianCoordinates &coordinates);
    ~CartesianCoordinates();

    // properties
    void resize(size_t size);
    size_t size() const;
    bool isEmpty() const;
    Matrix toMatrix() const;

    // coordinates
    void setPosition(size_t index, const Point3 &position);
    void setPosition(size_t index, Real x, Real y, Real z);
    Point3 position(size_t index) const;
    void setValue(size_t row, size_t column, Real value);
    Real value(size_t row, size_t column) const;
    void append(const Point3 &position);
    void append(Real x, Real y, Real z);
    void insert(size_t index, const Point3 &position);
    void insert(size_t index, Real x, Real y, Real z);
    void remove(size_t index);

    // geometry
    Real distance(size_t i, size_t j) const;
    Real angle(size_t i, size_t j, size_t k) const;
    Real angleRadians(size_t i, size_t j, size_t k) const;
    Real torsionAngle(size_t i, size_t j, size_t k, size_t l) const;
    Real torsionAngleRadians(size_t i, size_t j, size_t k, size_t l) const;
    Real wilsonAngle(size_t i, size_t j, size_t k, size_t l) const;
    Real wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l) const;
    Point3 center() const;
    Point3 weightedCenter(const std::vector<Real> &weights) const;
    void moveBy(const Vector3 &vector);
    void moveBy(Real x, Real y, Real z);
    Matrix distanceMatrix() const;

    // derivatives
    boost::array<Vector3, 2> distanceGradient(size_t i, size_t j) const;
    boost::array<Vector3, 3> angleGradient(size_t i, size_t j, size_t k) const;
    boost::array<Vector3, 3> angleGradientRadians(size_t i, size_t j, size_t k) const;
    boost::array<Vector3, 4> torsionAngleGradient(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> torsionAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> wilsonAngleGradient(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> wilsonAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const;

    // math
    CartesianCoordinates add(const CartesianCoordinates &coordinates) const;
    CartesianCoordinates subtract(const CartesianCoordinates &coordinates) const;
    Eigen::Matrix<Real, 3, 3> multiply(const CartesianCoordinates *coordinates) const;

    // operators
    CartesianCoordinates operator+(const CartesianCoordinates &coordinates) const;
    CartesianCoordinates operator-(const CartesianCoordinates &coordinates) const;
    CartesianCoordinates& operator=(const CartesianCoordinates &coordinates);

private:
    std::vector<Point3> m_coordinates;
};

} // end chemkit namespace

#endif // CHEMKIT_CARTESIANCOORDINATES_H
