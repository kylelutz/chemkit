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

template<typename T>
T function_pointer_cast(void *address)
{
    union FunctionPointer {
        void *address;
        T function;
    };

    FunctionPointer pointer;
    pointer.address = address;

    return pointer.function;
}

void setupGlFunctions()
{
    const QGLContext *context = QGLContext::currentContext();
    if(!context){
        return;
    }

    if(!glGenBuffers){
        glGenBuffers = function_pointer_cast<PFNGLGENBUFFERSPROC>(context->getProcAddress("glGenBuffersARB"));

        if(!glGenBuffers){
            return;
        }
    }

    if(!glBindBuffer){
        glBindBuffer = function_pointer_cast<PFNGLBINDBUFFERPROC>(context->getProcAddress("glBindBufferARB"));

        if(!glBindBuffer){
            return;
        }
    }

    if(!glBufferData){
        glBufferData = function_pointer_cast<PFNGLBUFFERDATAPROC>(context->getProcAddress("glBufferDataARB"));

        if(!glBufferData){
            return;
        }
    }

    if(!glBufferSubData){
        glBufferSubData = function_pointer_cast<PFNGLBUFFERSUBDATAPROC>(context->getProcAddress("glBufferSubDataARB"));

        if(!glBufferSubData){
            return;
        }
    }

    if(!glDeleteBuffers){
        glDeleteBuffers = function_pointer_cast<PFNGLDELETEBUFFERSPROC>(context->getProcAddress("glDeleteBuffersARB"));

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
    QVector<Vector3f> normals;
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
void GraphicsVertexBuffer::setNormals(const QVector<Vector3f> &normals)
{
    d->normals = normals;
}

/// Returns a list containing the vertex normals in the buffer.
QVector<Vector3f> GraphicsVertexBuffer::normals() const
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
                 d->verticies.size() * sizeof(Point3f) + d->normals.size() * sizeof(Vector3f),
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
                    d->normals.size() * sizeof(Vector3f),
                    d->normals.data());

    d->readyToDraw = true;
}

bool GraphicsVertexBuffer::readyToDraw() const
{
    return d->readyToDraw;
}

} // end chemkit namespace
