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

#include "scalarfield.h"

namespace chemkit {

// === ScalarFieldPrivate ================================================== //
class ScalarFieldPrivate
{
public:
    Point3 origin;
    std::vector<int> dimensions;
    std::vector<Real> lengths;
    std::vector<Real> data;
};

// === ScalarField ========================================================= //
/// \class ScalarField scalarfield.h chemkit/scalarfield.h
/// \ingroup chemkit
/// \brief The ScalarField class contains a three-dimensional grid of
///        scalar values.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty scalar field.
ScalarField::ScalarField()
    : d(new ScalarFieldPrivate)
{
    d->origin = Point3(0, 0, 0);
    d->dimensions = std::vector<int>(3, 0);
    d->lengths = std::vector<Real>(3, 0);
}

/// Creates a new scalar field.
ScalarField::ScalarField(const std::vector<int> &dimensions, const std::vector<Real> &cellLengths, const std::vector<Real> &data)
    : d(new ScalarFieldPrivate)
{
    d->origin = Point3(0, 0, 0);
    d->dimensions = dimensions;
    d->lengths = cellLengths;
    d->data = data;
}

/// Destroys the scalar field.
ScalarField::~ScalarField()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the width of the scalar field.
int ScalarField::width() const
{
    return d->dimensions[0];
}

/// Returns the height of the scalar field.
int ScalarField::height() const
{
    return d->dimensions[1];
}

/// Returns the depth of the scalar field.
int ScalarField::depth() const
{
    return d->dimensions[2];
}

/// Returns the size of the scalar field.
int ScalarField::size() const
{
    return width() * height() * depth();
}

/// Returns the dimensions of the scalar field.
std::vector<int> ScalarField::dimensions() const
{
    return d->dimensions;
}

/// Returns the width of a single cell in the grid.
Real ScalarField::cellWidth() const
{
    return d->lengths[0];
}

/// Returns the height of a single cell in the grid.
Real ScalarField::cellHeight() const
{
    return d->lengths[1];
}

/// Returns the depth of a single cell in the grid.
Real ScalarField::cellDepth() const
{
    return d->lengths[2];
}

/// Returns the dimensions of a single cell in the grid.
std::vector<Real> ScalarField::cellDimensions() const
{
    return d->lengths;
}

/// Sets the origin of the scalar field to \p origin.
void ScalarField::setOrigin(const Point3 &origin)
{
    d->origin = origin;
}

/// Returns the origin of the scalar field.
Point3 ScalarField::origin() const
{
    return d->origin;
}

/// Returns the data values for the scalar field.
std::vector<Real> ScalarField::data() const
{
    return d->data;
}

// --- Values -------------------------------------------------------------- //
/// Sets the value at (\p i, \p j, \p k) to \p value.
void ScalarField::setValue(int i, int j, int k, Real value)
{
    unsigned int index = i * d->dimensions[1] * d->dimensions[2] + j * d->dimensions[2] + k;

    if(index < d->data.size()){
        d->data[index] = value;
    }
}

/// Returns the the value at (\p i, \p j, \p k).
Real ScalarField::value(int i, int j, int k) const
{
    unsigned int index = i * d->dimensions[1] * d->dimensions[2] + j * d->dimensions[2] + k;

    if(index >= d->data.size()){
        return 0;
    }

    return d->data[index];
}

/// Returns the value at the position relative to the origin.
Real ScalarField::value(const Point3 &position) const
{
    Real x = position.x();
    Real y = position.y();
    Real z = position.z();

    int i = floor(x / d->lengths[0]);
    Real xd = x / d->lengths[0] - i;
    int j = floor(y / d->lengths[1]);
    Real yd = y / d->lengths[1] - j;
    int k = floor(z / d->lengths[2]);
    Real zd = z / d->lengths[2] - k;

    Real i1 = value(i + 0, j + 0, k + 0) * (1 - zd) + value(i + 0, j + 0, k + 1) * zd;
    Real i2 = value(i + 0, j + 1, k + 0) * (1 - zd) + value(i + 0, j + 1, k + 1) * zd;
    Real j1 = value(i + 1, j + 0, k + 0) * (1 - zd) + value(i + 1, j + 0, k + 1) * zd;
    Real j2 = value(i + 1, j + 1, k + 0) * (1 - zd) + value(i + 1, j + 1, k + 1) * zd;

    Real w1 = i1 * (1 - yd) + i2 * yd;
    Real w2 = j1 * (1 - yd) + j2 * yd;

    return w1 * (1 - xd) + w2 * xd;
}

/// Returns the position at (\p i, \p j, \p k).
Point3 ScalarField::position(int i, int j, int k) const
{
    return Point3(i * d->lengths[0],
                  j * d->lengths[1],
                  k * d->lengths[2]);
}

/// Returns the gradient at (\p i, \p j, \p k).
Vector3 ScalarField::gradient(int i, int j, int k) const
{
    return gradient(position(i, j, k));
}

/// Returns the gradient at the position relative to the origin.
Vector3 ScalarField::gradient(const Point3 &position) const
{
    Real h = 1.0e-4;

    return Vector3((value(position + Vector3(-h, 0, 0)) - value(position + Vector3(h, 0, 0))) / (2.0 * h),
                   (value(position + Vector3(0, -h, 0)) - value(position + Vector3(0, h, 0))) / (2.0 * h),
                   (value(position + Vector3(0, 0, -h)) - value(position + Vector3(0, 0, h))) / (2.0 * h));
}

} // end chemkit namespace
