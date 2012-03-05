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

#include "graphicssphere.h"

#include <chemkit/foreach.h>
#include <chemkit/geometry.h>

#include "graphicsray.h"
#include "graphicsvertexbuffer.h"

namespace chemkit {

namespace {

// (1 + sqrt(5)) / 2
const float GoldenRatio = 1.61803399f;

const float s = sqrt(1 + GoldenRatio*GoldenRatio);
const float t = GoldenRatio / s;
const float one = 1 / s;

const float IcosahedronVerticies[] = {
    t, one, 0,
    -t, one, 0,
    t, -one, 0,
    -t, -one, 0,
    one, 0, t,
    one, 0, -t,
    -one, 0, t,
    -one, 0, -t,
    0, t, one,
    0, -t, one,
    0, t, -one,
    0, -t, -one
};

const int IcosahedronVertexCount = 12;

const unsigned int IcosahedronIndicies[] = {
    0, 8, 4, 1, 10, 7, 2, 9, 11, 7, 3, 1, 0, 5, 10,
    3, 9, 6, 3, 11, 9, 8, 6, 4, 2, 4, 9, 3, 7, 11,
    4, 2, 0, 9, 4, 6, 2, 11, 5, 0, 10, 8, 5, 0, 2,
    10, 5, 7, 1, 6, 8, 1, 8, 10, 6, 1, 3, 11, 7, 5
};

const int IcosahedronIndexCount = 60;

} // end anonymous namespace

// === GraphicsSphere ====================================================== //
/// \class GraphicsSphere graphicssphere.h chemkit/graphicssphere.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsSphere class represents a sphere.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new sphere with radius of \c 0.
GraphicsSphere::GraphicsSphere()
{
    m_radius = 0;
}

/// Creates a new sphere with \p radius.
GraphicsSphere::GraphicsSphere(float radius)
{
    m_radius = radius;
}

/// Destroys the sphere object.
GraphicsSphere::~GraphicsSphere()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the radius of the sphere to \p radius.
void GraphicsSphere::setRadius(const float radius)
{
    m_radius = radius;
}

/// Returns the radius of the sphere.
float GraphicsSphere::radius() const
{
    return m_radius;
}

// --- Intersection -------------------------------------------------------- //
bool GraphicsSphere::intersects(const GraphicsRay &ray, float *distance) const
{
    Q_UNUSED(ray);
    Q_UNUSED(distance);

    return false;
}

// --- Tesselation --------------------------------------------------------- //
GraphicsVertexBuffer* GraphicsSphere::tesselate(int subdivisions) const
{
    // setup initial verticies
    QVector<Point3f> verticies(IcosahedronVertexCount);
    for(int i = 0; i < IcosahedronVertexCount; i++){
        const float *v = &IcosahedronVerticies[i*3];
        Point3f point(v[0], v[1], v[2]);
        point *= radius() / point.norm();
        verticies[i] = point;
    }

    // set up initial triangle indices
    QVector<unsigned short> indices(IcosahedronIndexCount);
    for(int i = 0; i < IcosahedronIndexCount; i++){
        indices[i] = IcosahedronIndicies[i];
    }

    // subdivide
    int subdivisionCount = 0;
    while(subdivisionCount < subdivisions){
        QVector<unsigned short> subdivisionIndicies;

        for(int i = 0; i < indices.size() / 3; i++){
            // indices of current triangle
            int i0 = indices[i*3];
            int i1 = indices[i*3+1];
            int i2 = indices[i*3+2];

            // verticies of the current triangle
            Point3f v0 = verticies[i0];
            Point3f v1 = verticies[i1];
            Point3f v2 = verticies[i2];

            // add three new verticies
            Point3f v01 = chemkit::geometry::midpoint(v0.cast<Real>(), v1.cast<Real>()).cast<float>();
            Point3f v12 = chemkit::geometry::midpoint(v1.cast<Real>(), v2.cast<Real>()).cast<float>();
            Point3f v20 = chemkit::geometry::midpoint(v2.cast<Real>(), v0.cast<Real>()).cast<float>();

            // scale points to lie on the sphere
            v01 *= radius() / v01.norm();
            v12 *= radius() / v12.norm();
            v20 *= radius() / v20.norm();

            // add verticies and record their indices
            int i01 = verticies.size();
            verticies.append(v01);
            int i12 = verticies.size();
            verticies.append(v12);
            int i20 = verticies.size();
            verticies.append(v20);

            // add triangle (i0, i01, i20)
            subdivisionIndicies.append(i0);
            subdivisionIndicies.append(i01);
            subdivisionIndicies.append(i20);

            // add triangle (i01, i12, i20)
            subdivisionIndicies.append(i01);
            subdivisionIndicies.append(i12);
            subdivisionIndicies.append(i20);

            // add triangle (i01, i1, i12)
            subdivisionIndicies.append(i01);
            subdivisionIndicies.append(i1);
            subdivisionIndicies.append(i12);

            // add triangle (i20, i12, i2)
            subdivisionIndicies.append(i20);
            subdivisionIndicies.append(i12);
            subdivisionIndicies.append(i2);
        }

        indices = subdivisionIndicies;

        subdivisionCount++;
    }

    // calculate vertex normals
    QVector<Vector3f> normals;
    foreach(const Point3f &vertex, verticies){
        normals.append(vertex.normalized());
    }

    // create vertex buffer
    GraphicsVertexBuffer *buffer = new GraphicsVertexBuffer;

    buffer->setVerticies(verticies);
    buffer->setNormals(normals);
    buffer->setIndicies(indices);

    return buffer;
}

} // end chemkit namespace
