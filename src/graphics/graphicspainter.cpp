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

#include "graphicspainter.h"

#include "graphicssphere.h"
#include "graphicscylinder.h"
#include "graphicsmaterial.h"
#include "graphicsvertexbuffer.h"

#include <chemkit/geometry.h>
#include <chemkit/constants.h>
#include <chemkit/quaternion.h>

#if defined(Q_WS_MAC)
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

namespace {

void nurbsErrorCallback(GLenum errorCode)
{
    qDebug() << "GraphicsPainter: Nurbs Error:" << (const char *) gluErrorString(errorCode);
}

} // end anonymous namespace

namespace chemkit {

// === GraphicsPainterPrivate ============================================== //
class GraphicsPainterPrivate
{
public:
    QColor drawColor;
};

// === GraphicsPainter ===================================================== //
/// \class GraphicsPainter graphicspainter.h chemkit/graphicspainter.h
/// \ingroup chemkit-graphics
/// \internal
/// \brief The GraphicsPainter class implements drawing methods.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics painter object.
GraphicsPainter::GraphicsPainter()
    : d(new GraphicsPainterPrivate)
{
}

/// Destroys the graphics painter object.
GraphicsPainter::~GraphicsPainter()
{
    delete d;
}

// --- Drawing ------------------------------------------------------------- //
void GraphicsPainter::draw(const GraphicsVertexBuffer *buffer, GraphicsPainter::PrimitiveType type)
{
    GLenum mode;

    switch(type){
    case Triangles:
        mode = GL_TRIANGLES;
        break;
    case TriangleStrip:
        mode = GL_TRIANGLE_STRIP;
        break;
    case TriangleFan:
        mode = GL_TRIANGLE_FAN;
        break;
    case Lines:
        mode = GL_LINES;
        break;
    case LineStrip:
        mode = GL_LINE_STRIP;
        break;
    case LineLoop:
        mode = GL_LINE_LOOP;
        break;
    case Points:
        mode = GL_POINTS;
        break;
    default:
        mode = GL_TRIANGLES;
        break;
    }

    buffer->draw(mode);
}

void GraphicsPainter::drawSphere(float radius)
{
    GraphicsSphere sphere(radius);
    GraphicsVertexBuffer *buffer = sphere.tesselate();

    draw(buffer);

    delete buffer;
}

void GraphicsPainter::drawSphere(const Point3f &center, float radius)
{
    glPushMatrix();
    glTranslated(center.x(), center.y(), center.z());

    drawSphere(radius);

    glPopMatrix();
}

void GraphicsPainter::drawCylinder(float radius, float length)
{
    GraphicsCylinder cylinder(radius, length);
    GraphicsVertexBuffer *buffer = cylinder.tesselate(12, 10);

    draw(buffer);

    delete buffer;
}

void GraphicsPainter::drawCylinder(const Point3f &a, const Point3f &b, float radius)
{
    glPushMatrix();

    glTranslatef(a.x(), a.y(), a.z());

    Vector3f vector = (a - b).normalized();
    Vector3f axis = vector.cross(-Vector3f::UnitZ()).normalized();
    float angle = chemkit::geometry::angle(vector.cast<Real>(), -Vector3f::UnitZ().cast<Real>());
    glRotatef(-angle, axis.x(), axis.y(), axis.z());

    float length = chemkit::geometry::distance(a.cast<Real>(), b.cast<Real>());

    drawCylinder(radius, length);

    glPopMatrix();
}

void GraphicsPainter::drawCircle(float radius)
{
    Q_UNUSED(radius);
}

void GraphicsPainter::drawCircle(const Point3f &center, float radius, const Vector3f &normal)
{
    Vector3f right(normal.y(), -normal.x(), 0); // vector orthogonal to normal

    glBegin(GL_TRIANGLE_FAN);

    glNormal3f(normal.x(), normal.y(), normal.z());
    glVertex3f(center.x(), center.y(), center.z());

    for(float angle = 0; angle <= 360; angle += 10){
        Point3f point = center + (chemkit::geometry::rotate(right, normal, angle).normalized() * radius);

        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(point.x(), point.y(), point.z());
    }

    glEnd();
}

void GraphicsPainter::drawTriangle(const Point3f &a, const Point3f &b, const Point3f &c)
{
    // verticies
    QVector<Point3f> verticies(3);
    verticies[0] = a;
    verticies[1] = b;
    verticies[2] = c;

    // normals
    Vector3f normal = (b - a).cross(c - b).normalized();
    QVector<Vector3f> normals(3);
    normals.fill(normal);

    // indices
    QVector<unsigned short> indices;
    indices << 0 << 1 << 2;

    // setup buffer
    GraphicsVertexBuffer buffer;
    buffer.setVerticies(verticies);
    buffer.setNormals(normals);
    buffer.setIndicies(indices);

    // draw buffer
    draw(&buffer);
}

void GraphicsPainter::drawRectangle(const Point3f &a, const Point3f &b, const Point3f &c, const Point3f &d)
{
    drawTriangle(a, b, c);
    drawTriangle(a, c, d);
}

void GraphicsPainter::drawSpline(const QList<Point3f> &points, float radius, int order)
{
    if(points.size() < 2){
        return;
    }

    // calculate axis and up vectors (needed to calculate control points)
    QVector<Vector3f> axisVectors(points.size());
    axisVectors[0] = (points[1] - points[0]).normalized();

    QVector<Vector3f> upVectors(points.size());
    if(points.size() > 2)
        upVectors[0] = (points[1] - points[0]).cross(points[2] - points[1]).normalized();
    else
        upVectors[0] = Vector3f::UnitZ();

    for(int i = 1; i < points.size(); i++){
        Vector3f axis = points[i] - points[i-1];

        if(i != (points.size() - 1)){
            float angle = chemkit::geometry::angle(axis.cast<Real>(), (points[i+1] - points[i]).cast<Real>());
            Vector3f rotationAxis = (points[i] - points[i-1]).cross(points[i+1] - points[i]).normalized();
            axis = chemkit::geometry::rotate(axis, rotationAxis, angle / 2.0f);
        }

        axisVectors[i] = axis.normalized();

        Vector3f rotationAxis = axisVectors[i-1].cross(axisVectors[i]);
        float angle = chemkit::geometry::angle(axis.cast<Real>(), axisVectors[i-1].cast<Real>());
        Vector3f up = chemkit::geometry::rotate(upVectors[i-1], rotationAxis, angle);
        up.normalize();

        upVectors[i] = up;
    }

    // calculate control points
    QVector<Point3f> controlPoints(points.size() * 9);

    for(int i = 0; i < points.size(); i++){
        const Point3f &point = points.at(i);

        Vector3f upVector = upVectors[i];
        Vector3f rightVector = upVector.cross(axisVectors[i]).normalized();

        // 8 points around a square surrounding point

        // right
        Point3f right = point + (rightVector * radius);
        controlPoints[i*9+0] = right;

        // bottom right
        Point3f bottomRight = right + (upVector * -radius);
        controlPoints[i*9+1] = bottomRight;

        // bottom
        Point3f bottom = point + (upVector * -radius);
        controlPoints[i*9+2] = bottom;

        // bottom left
        Point3f bottomLeft = bottom + (rightVector * -radius);
        controlPoints[i*9+3] = bottomLeft;

        // left
        Point3f left = point + (rightVector * -radius);
        controlPoints[i*9+4] = left;

        // top left
        Point3f topLeft = left + (upVector * radius);
        controlPoints[i*9+5] = topLeft;

        // top
        Point3f top = point + (upVector * radius);
        controlPoints[i*9+6] = top;

        // top right
        Point3f topRight = top + (rightVector * radius);

        controlPoints[i*9+7] = topRight;

        // right (again)
        controlPoints[i*9+8] = right;
    }

    // build knot vector
    QVector<float> uKnots(points.size() + order);
    for(int i = 0; i < uKnots.size(); i++){
        if(i < order)
            uKnots[i] = 0;
        else if(i > (points.size()) - 1)
            uKnots[i] = (points.size() + order) - 2 * order + 1;
        else
            uKnots[i] = (i - order + 1);
    }

    QVector<float> vKnots(12);
    vKnots[0] = 0;
    vKnots[1] = 0;
    vKnots[2] = 0;
    vKnots[3] = chemkit::constants::Pi / 2.0;
    vKnots[4] = chemkit::constants::Pi / 2.0;
    vKnots[5] = chemkit::constants::Pi;
    vKnots[6] = chemkit::constants::Pi;
    vKnots[7] = 3.0 * chemkit::constants::Pi / 2.0;
    vKnots[8] = 3.0 * chemkit::constants::Pi / 2.0;
    vKnots[9] = 2.0 * chemkit::constants::Pi;
    vKnots[10] = 2.0 * chemkit::constants::Pi;
    vKnots[11] = 2.0 * chemkit::constants::Pi;

    drawNurbsSurface(controlPoints, uKnots, vKnots, order, 3);

    drawCircle(points.first(), radius * 1.08, -axisVectors.first());
    drawCircle(points.last(), radius * 1.08, axisVectors.last());
}

void GraphicsPainter::drawNurbsSurface(const QVector<Point3f> &controlPoints, const QVector<float> &uKnots, const QVector<float> &vKnots, int uOrder, int vOrder)
{
    // nurbs control points
    QVector<float> points(controlPoints.size() * 3);
    for(int i = 0; i < controlPoints.size(); i++){
        points[i*3+0] = controlPoints[i].x();
        points[i*3+1] = controlPoints[i].y();
        points[i*3+2] = controlPoints[i].z();
    }

    // nurbs knot vectors
    QVector<float> uKnotVector(uKnots.size());
    for(int i = 0; i < uKnots.size(); i++)
        uKnotVector[i] = uKnots[i];

    QVector<float> vKnotVector(vKnots.size());
    for(int i = 0; i < vKnots.size(); i++)
        vKnotVector[i] = vKnots[i];

    int uStride = 3 * (vKnotVector.size() - vOrder);
    int vStride = 3;

    glEnable(GL_AUTO_NORMAL);
    GLUnurbsObj *nurb = gluNewNurbsRenderer();
    gluNurbsProperty(nurb, GLU_CULLING, GL_TRUE); // only render visible parts of the surface
    gluNurbsProperty(nurb, GLU_V_STEP, 4);
    gluNurbsProperty(nurb, GLU_U_STEP, 10);
    gluNurbsProperty(nurb, GLU_SAMPLING_METHOD, GLU_DOMAIN_DISTANCE);

#ifdef Q_OS_WIN
    gluNurbsCallback(nurb, GLU_ERROR, (void (__stdcall *)()) nurbsErrorCallback);
#else
    gluNurbsCallback(nurb, GLU_ERROR, (void (*)()) nurbsErrorCallback);
#endif

    gluBeginSurface(nurb);

    gluNurbsSurface(nurb,
                    uKnotVector.size(),
                    uKnotVector.data(),
                    vKnotVector.size(),
                    vKnotVector.data(),
                    uStride,
                    vStride,
                    points.data(),
                    uOrder,
                    vOrder,
                    GL_MAP2_VERTEX_3);

    gluEndSurface(nurb);
    gluDeleteNurbsRenderer(nurb);
    glDisable(GL_AUTO_NORMAL);
}

void GraphicsPainter::drawText(const QString &text, const QFont &font)
{
    Q_UNUSED(text);
    Q_UNUSED(font);
}

void GraphicsPainter::setColor(const QColor &color)
{
    d->drawColor = color;
    glColor4f(color.redF(), color.greenF(), color.blueF(), color.alphaF());
}

void GraphicsPainter::setMaterial(const GraphicsMaterial *material)
{
    glMateriali(GL_FRONT, GL_SHININESS, material->shininess());

    float specular[] = {static_cast<float>(material->specularColor().redF()),
                        static_cast<float>(material->specularColor().greenF()),
                        static_cast<float>(material->specularColor().blueF()),
                        static_cast<float>(material->specularColor().alphaF())};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
}

} // end chemkit namespace
