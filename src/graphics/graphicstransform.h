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
