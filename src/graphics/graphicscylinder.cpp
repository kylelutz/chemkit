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

#include "graphicscylinder.h"

#include "graphicsray.h"
#include "graphicsvertexbuffer.h"

namespace chemkit {

// === GraphicsCylinder ==================================================== //
/// \class GraphicsCylinder graphicscylinder.h chemkit/graphicscylinder.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsCylinder class represents a cylinder.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new cylinder with radius of \c 0 and length of \c 0.
GraphicsCylinder::GraphicsCylinder()
{
    m_radius = 0;
    m_length = 0;
}

/// Create a new cylinder object with \p radius and \p length.
GraphicsCylinder::GraphicsCylinder(float radius, float length)
{
    m_radius = radius;
    m_length = length;
}

/// Destroys the cylinder object.
GraphicsCylinder::~GraphicsCylinder()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the radius of the cylinder to \p radius.
void GraphicsCylinder::setRadius(float radius)
{
    m_radius = radius;
}

/// Returns the radius of the cylinder.
float GraphicsCylinder::radius() const
{
    return m_radius;
}

/// Sets the length of the cylinder to \p length.
void GraphicsCylinder::setLength(float length)
{
    m_length = length;
}

/// Returns the length of the cylinder.
float GraphicsCylinder::length() const
{
    return m_length;
}

// --- Intersection -------------------------------------------------------- //
bool GraphicsCylinder::intersects(const GraphicsRay &ray, float *distance) const
{
    Q_UNUSED(ray);
    Q_UNUSED(distance);

    return false;
}

// --- Tesselation --------------------------------------------------------- //
GraphicsVertexBuffer* GraphicsCylinder::tesselate(int slices, int stacks) const
{
    // slices must be at least 3
    slices = qMax(3, slices);

    // stacks must be at least 1
    stacks = qMax(1, stacks);

    float twoPi = chemkit::constants::Pi * 2;
    float sliceAngle = twoPi / slices;

    QVector<Point3f> verticies;

    for(int i = 0; i < slices; i++){
        float angle = i * sliceAngle;

        float x = cos(angle) * m_radius;
        float y = sin(angle) * m_radius;

        verticies.append(Point3f(x, y, 0));
    }

    Float stackHeight = m_length / stacks;

    for(int i = 1; i < stacks + 1; i++){
        for(int j = 0; j < slices; j++){
            Point3f lowerPoint = verticies[(i-1)*slices + j];

            verticies.append(lowerPoint.movedBy(0, 0, stackHeight));
        }
    }

    QVector<Vector3f> normals;
    foreach(const Point3f &point, verticies){
        Vector3f normal(point.x(), point.y(), 0);
        normal.normalize();
        normals.append(normal);
    }

    QVector<unsigned short> indicies;
    for(int i = 0; i < stacks; i++){
        for(int j = 0; j < slices; j++){
            // triangle indices
            unsigned short i0, i1, i2;

            // first triangle
            i0 = i * slices + j;
            i1 = i * slices + ((j + 1) % slices);
            i2 = (i+1) * slices + j;

            indicies.append(i0);
            indicies.append(i1);
            indicies.append(i2);

            // second triangle
            i0 = i1;
            i1 = (i+1)*slices + ((j + 1) % slices);

            indicies.append(i0);
            indicies.append(i1);
            indicies.append(i2);
        }
    }

    GraphicsVertexBuffer *buffer = new GraphicsVertexBuffer;

    buffer->setVerticies(verticies);
    buffer->setNormals(normals);
    buffer->setIndicies(indicies);

    return buffer;
}

} // end chemkit namespace
