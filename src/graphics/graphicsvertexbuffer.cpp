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

#include <chemkit/foreach.h>

#if !defined(Q_WS_MAC)
 #if defined(Q_WS_WIN)
   #include "../3rdparty/khronos/GL/glext.h"
 #else
  #include <GL/glext.h>
 #endif
#endif

#if !defined(Q_WS_MAC)
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
#endif // !Q_WS_MAC

namespace chemkit {

// === GraphicsVertexBufferPrivate ========================================= //
class GraphicsVertexBufferPrivate
{
public:
    bool readyToDraw;
    GLuint vertexBuffer;
    GLuint indexBuffer;
    QVector<Point3f> vertices;
    QVector<Vector3f> normals;
    QVector<unsigned short> indices;
    QVector<unsigned char> colors;
};

// === GraphicsVertexBuffer ================================================ //
/// \class GraphicsVertexBuffer graphicsvertexbuffer.h chemkit/graphicsvertexbuffer.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsVertexBuffer class represents a vertex buffer
///        object.
///
/// Vertex buffers contain vertex positions and optionally may also
/// contain data for normals, indices, and colors.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new, empty vertex buffer object.
GraphicsVertexBuffer::GraphicsVertexBuffer()
    : d(new GraphicsVertexBufferPrivate)
{
    d->readyToDraw = false;

#if !defined(Q_WS_MAC)
    setupGlFunctions();
#endif

    glGenBuffers(1, &d->vertexBuffer);
}

/// Create a new vertex buffer object and fill it with \p vertices.
GraphicsVertexBuffer::GraphicsVertexBuffer(const QVector<Point3f> &vertices)
    : d(new GraphicsVertexBufferPrivate)
{
    d->readyToDraw = false;
    d->vertices = vertices;

#if !defined(Q_WS_MAC)
    setupGlFunctions();
#endif
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
/// Returns the number of vertices in the buffer.
int GraphicsVertexBuffer::size() const
{
    return vertexCount();
}

/// Returns \c true if the vertex buffer is empty (i.e. size() == 0).
bool GraphicsVertexBuffer::isEmpty() const
{
    return size() == 0;
}

/// Removes all of the vertices and all associated data from the
/// buffer.
void GraphicsVertexBuffer::clear()
{
    d->vertices.clear();
    d->normals.clear();
    d->indices.clear();
    d->colors.clear();

    d->readyToDraw = false;
}

// --- Vertices ----------------------------------------------------------- //
/// Sets the vertices to \p vertices.
void GraphicsVertexBuffer::setVertices(const QVector<Point3f> &vertices)
{
    d->vertices = vertices;
}

/// Returns the vertices contained in the vertex buffer.
QVector<Point3f> GraphicsVertexBuffer::vertices() const
{
    return d->vertices;
}

/// Returns the number of vertices in the buffer.
int GraphicsVertexBuffer::vertexCount() const
{
    return d->vertices.size();
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

// --- Indices ------------------------------------------------------------- //
/// Sets the indices to \p indices.
void GraphicsVertexBuffer::setIndices(const QVector<unsigned short> &indices)
{
    d->indices = indices;
}

/// Returns the indices contained in the vertex buffer.
QVector<unsigned short> GraphicsVertexBuffer::indices() const
{
    return d->indices;
}

/// Returns the number of indices in the buffer.
int GraphicsVertexBuffer::indexCount() const
{
    return d->indices.size();
}

// --- Colors -------------------------------------------------------------- //
/// Sets the colors to \p colors.
void GraphicsVertexBuffer::setColors(const QVector<QColor> &colors)
{
    d->colors.clear();

    foreach(const QColor &color, colors){
        d->colors.append(color.red());
        d->colors.append(color.green());
        d->colors.append(color.blue());
        d->colors.append(color.alpha());
    }

    d->readyToDraw = false;
}

/// Returns the colors in the vertex buffer.
QVector<QColor> GraphicsVertexBuffer::colors() const
{
    QVector<QColor> colors;

    for(int i = 0; i < d->colors.size(); i += 4){
        colors.append(QColor::fromRgb(d->colors[i+0],
                                      d->colors[i+1],
                                      d->colors[i+2],
                                      d->colors[i+3]));
    }

    return colors;
}

/// Returns the number of colors in the vertex buffer.
int GraphicsVertexBuffer::colorCount() const
{
    return d->colors.size() / 4;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsVertexBuffer::draw(GLenum mode) const
{
    if(!d->readyToDraw){
        prepareToDraw();
    }

    glBindBuffer(GL_ARRAY_BUFFER, d->vertexBuffer);

    // setup vertices
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);

    // setup normals
    if(!d->normals.isEmpty()){
        size_t offset = d->vertices.size() * sizeof(Point3f);
        glNormalPointer(GL_FLOAT, 0, reinterpret_cast<void *>(offset));
        glEnableClientState(GL_NORMAL_ARRAY);
    }

    // setup colors
    if(!d->colors.isEmpty()){
        size_t offset = (d->vertices.size() * sizeof(Point3f)) +
                        (d->normals.size() * sizeof(Vector3f));
        glColorPointer(4, GL_UNSIGNED_BYTE, 0, reinterpret_cast<void *>(offset));
        glEnableClientState(GL_COLOR_ARRAY);
    }

    // draw
    if(!d->indices.isEmpty()){
        glDrawElements(mode, d->indices.size(), GL_UNSIGNED_SHORT, d->indices.data());
    }
    else{
        glDrawArrays(GL_POINTS, 0, d->vertices.size());
    }

    // cleanup state
    if(!d->colors.isEmpty()){
        glDisableClientState(GL_COLOR_ARRAY);
    }

    if(!d->normals.isEmpty()){
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
}

void GraphicsVertexBuffer::prepareToDraw() const
{
    glBindBuffer(GL_ARRAY_BUFFER, d->vertexBuffer);

    // allocate space
    glBufferData(GL_ARRAY_BUFFER,
                 d->vertices.size() * sizeof(Point3f) +
                 d->normals.size() * sizeof(Vector3f) +
                 d->colors.size() * sizeof(unsigned char),
                 0,
                 GL_STATIC_DRAW);

    // load vertices
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    d->vertices.size() * sizeof(Point3f),
                    d->vertices.data());

    // load normals
    glBufferSubData(GL_ARRAY_BUFFER,
                    d->vertices.size() * sizeof(Point3f),
                    d->normals.size() * sizeof(Vector3f),
                    d->normals.data());

    // load colors
    glBufferSubData(GL_ARRAY_BUFFER,
                    d->vertices.size() * sizeof(Point3f) +
                    d->normals.size() * sizeof(Vector3f),
                    d->colors.size() * sizeof(unsigned char),
                    d->colors.data());

    d->readyToDraw = true;
}

bool GraphicsVertexBuffer::readyToDraw() const
{
    return d->readyToDraw;
}

} // end chemkit namespace
