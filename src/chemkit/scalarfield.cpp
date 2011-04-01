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

#include "scalarfield.h"

namespace chemkit {

// === ScalarFieldPrivate ================================================== //
class ScalarFieldPrivate
{
    public:
        Point3 origin;
        std::vector<int> dimensions;
        std::vector<Float> lengths;
        std::vector<Float> data;
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
    d->dimensions = std::vector<int>(3, 0);
    d->lengths = std::vector<Float>(3, 0);
}

/// Creates a new scalar field.
ScalarField::ScalarField(const std::vector<int> &dimensions, const std::vector<Float> &cellLengths, const std::vector<Float> &data)
    : d(new ScalarFieldPrivate)
{
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
Float ScalarField::cellWidth() const
{
    return d->lengths[0];
}

/// Returns the height of a single cell in the grid.
Float ScalarField::cellHeight() const
{
    return d->lengths[1];
}

/// Returns the depth of a single cell in the grid.
Float ScalarField::cellDepth() const
{
    return d->lengths[2];
}

/// Returns the dimensions of a single cell in the grid.
std::vector<Float> ScalarField::cellDimensions() const
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
std::vector<Float> ScalarField::data() const
{
    return d->data;
}

// --- Values -------------------------------------------------------------- //
/// Sets the value at (\p i, \p j, \p k) to \p value.
void ScalarField::setValue(int i, int j, int k, Float value)
{
    unsigned int index = i * d->dimensions[1] * d->dimensions[2] + j * d->dimensions[2] + k;

    if(index < d->data.size()){
        d->data[index] = value;
    }
}

/// Returns the the value at (\p i, \p j, \p k).
Float ScalarField::value(int i, int j, int k) const
{
    unsigned int index = i * d->dimensions[1] * d->dimensions[2] + j * d->dimensions[2] + k;

    if(index >= d->data.size()){
        return 0;
    }

    return d->data[index];
}

/// Returns the value at the position relative to the origin.
Float ScalarField::value(const Point3 &position) const
{
    Float x = position.x();
    Float y = position.y();
    Float z = position.z();

    int i = floor(x / d->lengths[0]);
    Float xd = x / d->lengths[0] - i;
    int j = floor(y / d->lengths[1]);
    Float yd = y / d->lengths[1] - j;
    int k = floor(z / d->lengths[2]);
    Float zd = z / d->lengths[2] - k;

    Float i1 = value(i + 0, j + 0, k + 0) * (1 - zd) + value(i + 0, j + 0, k + 1) * zd;
    Float i2 = value(i + 0, j + 1, k + 0) * (1 - zd) + value(i + 0, j + 1, k + 1) * zd;
    Float j1 = value(i + 1, j + 0, k + 0) * (1 - zd) + value(i + 1, j + 0, k + 1) * zd;
    Float j2 = value(i + 1, j + 1, k + 0) * (1 - zd) + value(i + 1, j + 1, k + 1) * zd;

    Float w1 = i1 * (1 - yd) + i2 * yd;
    Float w2 = j1 * (1 - yd) + j2 * yd;

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
    Float h = 1.0e-4;

    return Vector3((value(position.movedBy(-h, 0, 0)) - value(position.movedBy(h, 0, 0))) / (2.0 * h),
                   (value(position.movedBy(0, -h, 0)) - value(position.movedBy(0, h, 0))) / (2.0 * h),
                   (value(position.movedBy(0, 0, -h)) - value(position.movedBy(0, 0, h))) / (2.0 * h));
}

} // end chemkit namespace
