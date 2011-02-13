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

#include "graphicssphere.h"

#include "graphicsray.h"
#include "graphicsvertexbuffer.h"

namespace chemkit {

namespace {

// (1 + sqrt(5)) / 2
const GraphicsFloat GoldenRatio = chemkit::GraphicsFloat(1.61803399);

const GraphicsFloat s = sqrt(1 + GoldenRatio*GoldenRatio);
const GraphicsFloat t = GoldenRatio / s;
const GraphicsFloat one = 1 / s;

const GraphicsFloat IcosahedronVerticies[] = {
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
GraphicsSphere::GraphicsSphere(GraphicsFloat radius)
{
    m_radius = radius;
}

/// Destroys the sphere object.
GraphicsSphere::~GraphicsSphere()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the radius of the sphere to \p radius.
void GraphicsSphere::setRadius(const GraphicsFloat radius)
{
    m_radius = radius;
}

/// Returns the radius of the sphere.
GraphicsFloat GraphicsSphere::radius() const
{
    return m_radius;
}

// --- Intersection -------------------------------------------------------- //
bool GraphicsSphere::intersects(const GraphicsRay &ray, GraphicsFloat *distance) const
{
    Q_UNUSED(ray);
    Q_UNUSED(distance);

    return false;
}

// --- Tesselation --------------------------------------------------------- //
GraphicsVertexBuffer* GraphicsSphere::tesselate(int subdivisions) const
{
    // setup initial verticies
    QVector<GraphicsPoint> verticies(IcosahedronVertexCount);
    for(int i = 0; i < IcosahedronVertexCount; i++){
        const GraphicsFloat *v = &IcosahedronVerticies[i*3];
        GraphicsPoint point(v[0], v[1], v[2]);
        point.scale(radius() / point.norm());
        verticies[i] = point;
    }

    // set up initial triangle indicies
    QVector<unsigned int> indicies(IcosahedronIndexCount);
    for(int i = 0; i < IcosahedronIndexCount; i++){
        indicies[i] = IcosahedronIndicies[i];
    }

    // subdivide
    int subdivisionCount = 0;
    while(subdivisionCount < subdivisions){
        QVector<unsigned int> subdivisionIndicies;

        for(int i = 0; i < indicies.size() / 3; i++){
            // indicies of current triangle
            int i0 = indicies[i*3];
            int i1 = indicies[i*3+1];
            int i2 = indicies[i*3+2];

            // verticies of the current triangle
            GraphicsPoint v0 = verticies[i0];
            GraphicsPoint v1 = verticies[i1];
            GraphicsPoint v2 = verticies[i2];

            // add three new verticies
            GraphicsPoint v01 = v0.midpoint(v1);
            GraphicsPoint v12 = v1.midpoint(v2);
            GraphicsPoint v20 = v2.midpoint(v0);

            // scale points to lie on the sphere
            v01.scale(radius() / v01.norm());
            v12.scale(radius() / v12.norm());
            v20.scale(radius() / v20.norm());

            // add verticies and record their indicies
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

        indicies = subdivisionIndicies;

        subdivisionCount++;
    }

    // calculate vertex normals
    QVector<GraphicsVector> normals;
    foreach(const GraphicsPoint &vertex, verticies){
        normals.append(GraphicsVector(vertex).normalized());
    }

    // create vertex buffer
    GraphicsVertexBuffer *buffer = new GraphicsVertexBuffer;

    buffer->setVerticies(verticies);
    buffer->setNormals(normals);
    buffer->setIndicies(indicies);

    return buffer;
}

} // end chemkit namespace
