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

#include "graphicsvertexbuffer.h"

#include <QtOpenGL>

#ifdef Q_OS_WIN32
#include "../3rdparty/khronos/GL/glext.h"
#else
#include <GL/glext.h>
#endif

namespace {

PFNGLGENBUFFERSPROC glGenBuffers = 0;
PFNGLBINDBUFFERPROC glBindBuffer = 0;
PFNGLBUFFERDATAPROC glBufferData = 0;
PFNGLBUFFERSUBDATAPROC glBufferSubData = 0;
PFNGLDELETEBUFFERSPROC glDeleteBuffers = 0;

void setupGlFunctions()
{
    const QGLContext *context = QGLContext::currentContext();
    if(!context){
        return;
    }

    if(!glGenBuffers){
        glGenBuffers = (PFNGLGENBUFFERSPROC) context->getProcAddress("glGenBuffersARB");

        if(!glGenBuffers){
            return;
        }
    }

    if(!glBindBuffer){
        glBindBuffer = (PFNGLBINDBUFFERPROC) context->getProcAddress("glBindBufferARB");

        if(!glBindBuffer){
            return;
        }
    }

    if(!glBufferData){
        glBufferData = (PFNGLBUFFERDATAPROC) context->getProcAddress("glBufferDataARB");

        if(!glBufferData){
            return;
        }
    }

    if(!glBufferSubData){
        glBufferSubData = (PFNGLBUFFERSUBDATAPROC) context->getProcAddress("glBufferSubDataARB");

        if(!glBufferSubData){
            return;
        }
    }

    if(!glDeleteBuffers){
        glDeleteBuffers = (PFNGLDELETEBUFFERSPROC) context->getProcAddress("glDeleteBuffersARB");

        if(!glDeleteBuffers){
            return;
        }
    }
}

} // end anonymous namespace

namespace chemkit {

// === GraphicsVertexBufferPrivate ========================================= //
class GraphicsVertexBufferPrivate
{
    public:
        bool readyToDraw;
        GLuint vertexBuffer;
        GLuint indexBuffer;
        QVector<Point3f> verticies;
        QVector<Vector3g> normals;
        QVector<unsigned short> indicies;
};

// === GraphicsVertexBuffer ================================================ //
/// \class GraphicsVertexBuffer graphicsvertexbuffer.h chemkit/graphicsvertexbuffer.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsVertexBuffer class represents a vertex buffer
///        object.
///
/// Vertex buffers contain vertex positions and optionally may also
/// contain data for normals, indicies, and colors.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new, empty vertex buffer object.
GraphicsVertexBuffer::GraphicsVertexBuffer()
    : d(new GraphicsVertexBufferPrivate)
{
    d->readyToDraw = false;

    setupGlFunctions();

    glGenBuffers(1, &d->vertexBuffer);
}

/// Create a new vertex buffer object and fill it with \p verticies.
GraphicsVertexBuffer::GraphicsVertexBuffer(const QVector<Point3f> &verticies)
    : d(new GraphicsVertexBufferPrivate)
{
    d->readyToDraw = false;
    d->verticies = verticies;

    setupGlFunctions();
    glGenBuffers(1, &d->vertexBuffer);
}

/// Destroys the graphics vertex buffer object.
GraphicsVertexBuffer::~GraphicsVertexBuffer()
{
    if(d->vertexBuffer){
        glDeleteBuffers(1, &d->vertexBuffer);
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of verticies in the buffer.
int GraphicsVertexBuffer::size() const
{
    return vertexCount();
}

/// Returns \c true if the vertex buffer is empty (i.e. size() == 0).
bool GraphicsVertexBuffer::isEmpty() const
{
    return size() == 0;
}

/// Removes all of the verticies and all associated data from the
/// buffer.
void GraphicsVertexBuffer::clear()
{
    d->verticies.clear();
    d->normals.clear();
    d->indicies.clear();

    d->readyToDraw = false;
}

// --- Verticies ----------------------------------------------------------- //
/// Sets the verticies to \p verticies.
void GraphicsVertexBuffer::setVerticies(const QVector<Point3f> &verticies)
{
    d->verticies = verticies;
}

/// Returns the verticies contained in the vertex buffer.
QVector<Point3f> GraphicsVertexBuffer::verticies() const
{
    return d->verticies;
}

/// Returns the number of verticies in the buffer.
int GraphicsVertexBuffer::vertexCount() const
{
    return d->verticies.size();
}

// --- Normals ------------------------------------------------------------- //
/// Sets the vertex normals to \p normals.
void GraphicsVertexBuffer::setNormals(const QVector<Vector3g> &normals)
{
    d->normals = normals;
}

/// Returns a list containing the vertex normals in the buffer.
QVector<Vector3g> GraphicsVertexBuffer::normals() const
{
    return d->normals;
}

/// Returns the number of vertex normals in the buffer.
int GraphicsVertexBuffer::normalCount() const
{
    return d->normals.size();
}

// --- Indicies ------------------------------------------------------------ //
/// Sets the indicies to \p indicies.
void GraphicsVertexBuffer::setIndicies(const QVector<unsigned short> &indicies)
{
    d->indicies = indicies;
}

/// Returns the indicies contained in the vertex buffer.
QVector<unsigned short> GraphicsVertexBuffer::indicies() const
{
    return d->indicies;
}

/// Returns the number of indicies in the buffer.
int GraphicsVertexBuffer::indexCount() const
{
    return d->indicies.size();
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsVertexBuffer::draw() const
{
    if(!d->readyToDraw){
        prepareToDraw();
    }

    glBindBuffer(GL_ARRAY_BUFFER, d->vertexBuffer);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glNormalPointer(GL_FLOAT, 0, reinterpret_cast<void *>(d->verticies.size() * sizeof(Point3f)));

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    if(!d->indicies.isEmpty()){
        glDrawElements(GL_TRIANGLES, d->indicies.size(), GL_UNSIGNED_SHORT, d->indicies.data());
    }
    else{
        glDrawArrays(GL_POINTS, 0, d->verticies.size());
    }
}

void GraphicsVertexBuffer::prepareToDraw() const
{
    glBindBuffer(GL_ARRAY_BUFFER, d->vertexBuffer);

    // allocate space
    glBufferData(GL_ARRAY_BUFFER,
                 d->verticies.size() * sizeof(Point3f) + d->normals.size() * sizeof(Vector3g),
                 0,
                 GL_STATIC_DRAW);

    // load verticies
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    d->verticies.size() * sizeof(Point3f),
                    d->verticies.data());

    // load normals
    glBufferSubData(GL_ARRAY_BUFFER,
                    d->verticies.size() * sizeof(Point3f),
                    d->normals.size() * sizeof(Vector3g),
                    d->normals.data());

    d->readyToDraw = true;
}

bool GraphicsVertexBuffer::readyToDraw() const
{
    return d->readyToDraw;
}

} // end chemkit namespace
