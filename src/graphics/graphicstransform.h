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

#ifndef CHEMKIT_GRAPHICSTRANSFORM_H
#define CHEMKIT_GRAPHICSTRANSFORM_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>
#include <chemkit/staticmatrix.h>

#include "graphicsray.h"

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsTransform
{
    public:
        // construction and destruction
        GraphicsTransform();
        GraphicsTransform(const GraphicsTransform &transform);
        GraphicsTransform(const StaticMatrix<float, 4, 4> &matrix);
        ~GraphicsTransform();

        // properties
        const float* data() const;

        // math
        void invert();
        GraphicsTransform inverted() const;
        GraphicsRay multiply(const GraphicsRay &ray) const;
        Point3f multiply(const Point3f &point) const;
        Vector3f multiply(const Vector3f &vector) const;
        StaticVector<float, 4> multiply(const StaticVector<float, 4> &vector);
        GraphicsTransform multiply(const GraphicsTransform &transform) const;
        Point3f inverseMultiply(const Point3f &point) const;
        Vector3f inverseMultiply(const Vector3f &vector) const;
        StaticVector<float, 4> inverseMultiply(const StaticVector<float, 4> &vector);

        // operators
        float operator()(int row, int column) const;
        float& operator()(int row, int column);
        GraphicsRay operator*(const GraphicsRay &ray) const;
        Point3f operator*(const Point3f &point) const;
        Vector3f operator*(const Vector3f &vector) const;
        GraphicsTransform operator*(const GraphicsTransform &transform) const;
        GraphicsTransform& operator*=(const GraphicsTransform &transform);
        GraphicsTransform& operator=(const GraphicsTransform &transform);
        CommaInitializer<float> operator<<(const float value);

        // static methods
        static GraphicsTransform identity();
        static GraphicsTransform translation(const Vector3f &vector);
        static GraphicsTransform rotation(const Vector3f &axis, float angle);
        static GraphicsTransform perspective(float angle, float aspectRatio, float nearDistance, float farDistance);
        static GraphicsTransform frustum(float left, float right, float top, float bottom, float nearDistance, float farDistance);
        static GraphicsTransform orthographic(float left, float right, float top, float bottom, float near, float far);

    private:
        StaticMatrix<float, 4, 4> *m_matrix;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSTRANSFORM_H
