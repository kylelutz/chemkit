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

#include "graphicstransform.h"

namespace chemkit {

// === GraphicsTransform =================================================== //
/// \class GraphicsTransform graphicstransform.h chemkit/graphicstransform.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsTransform class represents a transformation
///        matrix.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty graphics transform.
///
/// The tranformation returned is:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       0 & 0 & 0 & 0 \\
///       0 & 0 & 0 & 0 \\
///       0 & 0 & 0 & 0 \\
///       0 & 0 & 0 & 0 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform::GraphicsTransform()
{
    m_matrix = new StaticMatrix<float, 4, 4>();
}

/// Creates a new transform as a copy of \p transform.
GraphicsTransform::GraphicsTransform(const GraphicsTransform &transform)
{
    m_matrix = new StaticMatrix<float, 4, 4>(*transform.m_matrix);
}

/// Creates a new transform that contains \p matrix.
GraphicsTransform::GraphicsTransform(const StaticMatrix<float, 4, 4> &matrix)
{
    m_matrix = new StaticMatrix<float, 4, 4>(matrix);
}

/// Destroys the graphics transform.
GraphicsTransform::~GraphicsTransform()
{
    delete m_matrix;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the data for the transform.
///
/// Use the following code to load a GraphicsTransform into
/// OpenGL:
/// \code
/// glLoadMatrixf(transform.data());
/// \endcode
const float* GraphicsTransform::data() const
{
    return m_matrix->data();
}

// --- Math ---------------------------------------------------------------- //
/// Inverts the transform.
void GraphicsTransform::invert()
{
    m_matrix->invert();
}

/// Returns the inverted version of the transform.
GraphicsTransform GraphicsTransform::inverted() const
{
    GraphicsTransform transform = *this;
    transform.invert();
    return transform;
}

/// Multiplies \p ray by the transform.
GraphicsRay GraphicsTransform::multiply(const GraphicsRay &ray) const
{
    Point3g origin = multiply(ray.origin());
    Point3g direction = multiply(ray.direction());

    return GraphicsRay(origin, direction);
}

/// Multiplies \p point by the transform.
Point3g GraphicsTransform::multiply(const Point3g &point) const
{
    StaticVector<float, 4> vector4;
    vector4[0] = point.x();
    vector4[1] = point.y();
    vector4[2] = point.z();
    vector4[3] = 1;

    vector4 = m_matrix->multiply(vector4);

    return Point3g(vector4[0], vector4[1], vector4[2]);
}

/// Multiplies \p vector by the transform.
Vector3g GraphicsTransform::multiply(const Vector3g &vector) const
{
    StaticVector<float, 4> vector4;
    vector4[0] = vector.x();
    vector4[1] = vector.y();
    vector4[2] = vector.z();
    vector4[3] = 0;

    vector4 = m_matrix->multiply(vector4);

    return Vector3g(vector4[0], vector4[1], vector4[2]);
}

/// Multiplies \p transform by the transform.
GraphicsTransform GraphicsTransform::multiply(const GraphicsTransform &transform) const
{
    return m_matrix->multiply(*transform.m_matrix);
}

StaticVector<float, 4> GraphicsTransform::multiply(const StaticVector<float, 4> &vector)
{
    return m_matrix->multiply(vector);
}

/// Multiplies \p point by the inverse of the transform.
Point3g GraphicsTransform::inverseMultiply(const Point3g &point) const
{
    StaticVector<float, 4> vector4;
    vector4[0] = point.x();
    vector4[1] = point.y();
    vector4[2] = point.z();
    vector4[3] = 1;

    vector4 = m_matrix->inverted().multiply(vector4);

    return Point3g(vector4[0], vector4[1], vector4[2]);
}

/// Multiplies \p vector by the inverse of the transform.
Vector3g GraphicsTransform::inverseMultiply(const Vector3g &vector) const
{
    StaticVector<float, 4> vector4;
    vector4[0] = vector.x();
    vector4[1] = vector.y();
    vector4[2] = vector.z();
    vector4[3] = 0;

    vector4 = m_matrix->inverted().multiply(vector4);

    return Vector3g(vector4[0], vector4[1], vector4[2]);
}

StaticVector<float, 4> GraphicsTransform::inverseMultiply(const StaticVector<float, 4> &vector)
{
    return m_matrix->inverted().multiply(vector);
}

// --- Operators ----------------------------------------------------------- //
float GraphicsTransform::operator()(int row, int column) const
{
    return m_matrix->operator()(row, column);
}

float& GraphicsTransform::operator()(int row, int column)
{
    return m_matrix->operator()(row, column);
}

GraphicsRay GraphicsTransform::operator*(const GraphicsRay &ray) const
{
    return multiply(ray);
}

Point3g GraphicsTransform::operator*(const Point3g &point) const
{
    return multiply(point);
}

Vector3g GraphicsTransform::operator*(const Vector3g &vector) const
{
    return multiply(vector);
}

GraphicsTransform GraphicsTransform::operator*(const GraphicsTransform &transform) const
{
    return multiply(transform);
}

GraphicsTransform& GraphicsTransform::operator*=(const GraphicsTransform &transform)
{
    *m_matrix *= *transform.m_matrix;
    return *this;
}

GraphicsTransform& GraphicsTransform::operator=(const GraphicsTransform &transform)
{
    *m_matrix = *transform.m_matrix;
    return *this;
}

CommaInitializer<float> GraphicsTransform::operator<<(const float value)
{
    m_matrix->data()[0] = value;

    return CommaInitializer<float>(m_matrix->data(), 4, 4);
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the identity transform.
///
/// The transformation returned is the following:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       1 & 0 & 0 & 0 \\
///       0 & 1 & 0 & 0 \\
///       0 & 0 & 1 & 0 \\
///       0 & 0 & 0 & 1 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform GraphicsTransform::identity()
{
    return StaticMatrix<float, 4, 4>::identity();
}

/// Returns a transformation matrix that represents the translation by
/// \p vector.
///
/// The transformation returned is the following:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       1 & 0 & 0 & vector_{x} \\
///       0 & 1 & 0 & vector_{y} \\
///       0 & 0 & 1 & vector_{z} \\
///       0 & 0 & 0 & 1 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform GraphicsTransform::translation(const Vector3g &vector)
{
    GraphicsTransform transform = identity();

    transform(0, 3) = vector.x();
    transform(1, 3) = vector.y();
    transform(2, 3) = vector.z();

    return transform;
}

/// Returns a transform that represents a rotation by \p angle
/// degrees around \p axis.
GraphicsTransform GraphicsTransform::rotation(const Vector3g &axis, float angle)
{
    GraphicsTransform transform = identity();

    Vector3g v = axis.normalized();
    float c = cos(angle * chemkit::constants::DegreesToRadians);
    float s = sin(angle * chemkit::constants::DegreesToRadians);

    transform(0, 0) = v.x() * v.x() + (1 - v.x() * v.x()) * c;
    transform(0, 1) = v.x() * v.y() * (1 - c) - v.z() * s;
    transform(0, 2) = v.x() * v.z() * (1 - c) + v.y() * s;
    transform(1, 0) = v.x() * v.y() * (1 - c) + v.z() * s;
    transform(1, 1) = v.y() * v.y() + (1 - v.y() * v.y()) * c;
    transform(1, 2) = v.y() * v.z() * (1 - c) - v.x() * s;
    transform(2, 0) = v.x() * v.z() * (1 - c) - v.y() * s;
    transform(2, 1) = v.y() * v.z() * (1 - c) + v.x() * s;
    transform(2, 2) = v.z() * v.z() + (1 - v.z() * v.z()) * c;

    return transform;
}

/// Returns a perspective transform.
///
/// The transformation returned is the following:
/// \f[ f = cot(\frac{angle}{2}) \f]
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       \frac{f}{aspectRatio} & 0 & 0 & 0 \\
///       0 & f & 0 & 0 \\
///       0 & 0 & \frac{nearDistance+farDistance}{nearDistance-farDistance} & \frac{2 \cdot nearDistance \cdot farDistance}{nearDistance-farDistance} \\
///       0 & 0 & -1 & 0 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform GraphicsTransform::perspective(float angle, float aspectRatio, float nearDistance, float farDistance)
{
    GraphicsTransform transform;

    float f = 1.0 / tan(angle / 2.0);

    transform(0, 0) = f / aspectRatio;
    transform(1, 1) = f;
    transform(2, 2) = (nearDistance + farDistance) / (nearDistance - farDistance);
    transform(2, 3) = (2 * nearDistance * farDistance) / (nearDistance - farDistance);
    transform(3, 2) = -1;

    return transform;
}

/// Returns a frustum transform.
///
/// The transformation returned is the following:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       \frac{2 \cdot nearDistance}{right-left} & 0 & \frac{right+left}{right-left} & 0 \\
///       0 & \frac{2 \cdot nearDistance}{top-bottom} & \frac{top+bottom}{top-bottom} & 0 \\
///       0 & 0 & -\frac{farDistance+nearDistance}{farDistance-nearDistance} & -\frac{2 \cdot farDistance \cdot nearDistance}{farDistance-nearDistance} \\
///       0 & 0 & -1 & 0 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform GraphicsTransform::frustum(float left, float right, float top, float bottom, float nearDistance, float farDistance)
{
    GraphicsTransform transform;

    transform(0, 0) = (2 * nearDistance) / (right - left);
    transform(1, 1) = (2 * nearDistance) / (top - bottom);
    transform(2, 0) = (right + left) / (right - left);
    transform(2, 1) = (top + bottom) / (top - bottom);
    transform(2, 2) = -(farDistance + nearDistance) / (farDistance - nearDistance);
    transform(2, 3) = -(2 * farDistance * nearDistance) / (farDistance - nearDistance);
    transform(3, 2) = -1;

    return transform;
}

/// Returns a orthographic transform.
///
/// The transformation returned is the following:
/** \f[
///   \left[
///   {
///     \begin{array}{cccc}
///       \frac{2}{right-left} & 0 & 0 & -\frac{right+left}{right-left} \\
///       0 & \frac{2}{top-bottom} & 0 & -\frac{top+bottom}{top-bottom} \\
///       0 & 0 & -\frac{2}{farDistance-nearDistance} & -\frac{farDistance+nearDistance}{farDistance-nearDistance} \\
///       0 & 0 & 0 & 1 \\
///     \end{array}
///   }
///   \right]
/// \f]
**/
GraphicsTransform GraphicsTransform::orthographic(float left, float right, float top, float bottom, float nearDistance, float farDistance)
{
    GraphicsTransform transform;

    transform(0, 0) = 2.0 / (right - left);
    transform(0, 3) = -(right + left) / (right - left);
    transform(1, 1) = 2.0 / (top - bottom);
    transform(1, 3) = -(top + bottom) / (top - bottom);
    transform(2, 2) = -2.0 / (farDistance - nearDistance);
    transform(2, 3) = -(farDistance + nearDistance) / (farDistance - nearDistance);
    transform(3, 3) = 1;

    return transform;
}

} // end chemkit namespace
