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

#include "point3g.h"
#include "vector3g.h"
#include "graphicsray.h"

#include <chemkit/staticmatrix.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsTransform
{
    public:
        // construction and destruction
        GraphicsTransform();
        GraphicsTransform(const GraphicsTransform &transform);
        GraphicsTransform(const StaticMatrix<GraphicsFloat, 4, 4> &matrix);
        ~GraphicsTransform();

        // properties
        const GraphicsFloat* data() const;

        // math
        void invert();
        GraphicsTransform inverted() const;
        GraphicsRay multiply(const GraphicsRay &ray) const;
        Point3g multiply(const Point3g &point) const;
        Vector3g multiply(const Vector3g &vector) const;
        StaticVector<GraphicsFloat, 4> multiply(const StaticVector<GraphicsFloat, 4> &vector);
        GraphicsTransform multiply(const GraphicsTransform &transform) const;
        Point3g inverseMultiply(const Point3g &point) const;
        Vector3g inverseMultiply(const Vector3g &vector) const;
        StaticVector<GraphicsFloat, 4> inverseMultiply(const StaticVector<GraphicsFloat, 4> &vector);

        // operators
        GraphicsFloat operator()(int row, int column) const;
        GraphicsFloat& operator()(int row, int column);
        GraphicsRay operator*(const GraphicsRay &ray) const;
        Point3g operator*(const Point3g &point) const;
        Vector3g operator*(const Vector3g &vector) const;
        GraphicsTransform operator*(const GraphicsTransform &transform) const;
        GraphicsTransform& operator*=(const GraphicsTransform &transform);
        GraphicsTransform& operator=(const GraphicsTransform &transform);
        CommaInitializer<GraphicsFloat> operator<<(const GraphicsFloat value);

        // static methods
        static GraphicsTransform identity();
        static GraphicsTransform translation(const Vector3g &vector);
        static GraphicsTransform rotation(const Vector3g &axis, GraphicsFloat angle);
        static GraphicsTransform perspective(GraphicsFloat angle, GraphicsFloat aspectRatio, GraphicsFloat nearDistance, GraphicsFloat farDistance);
        static GraphicsTransform frustum(GraphicsFloat left, GraphicsFloat right, GraphicsFloat top, GraphicsFloat bottom, GraphicsFloat nearDistance, GraphicsFloat farDistance);
        static GraphicsTransform orthographic(GraphicsFloat left, GraphicsFloat right, GraphicsFloat top, GraphicsFloat bottom, GraphicsFloat near, GraphicsFloat far);

    private:
        StaticMatrix<GraphicsFloat, 4, 4> *m_matrix;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSTRANSFORM_H
