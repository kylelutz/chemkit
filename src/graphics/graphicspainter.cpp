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

#include "graphicspainter.h"

#include "graphicssphere.h"
#include "graphicscylinder.h"
#include "graphicsmaterial.h"
#include "graphicsquaternion.h"
#include "graphicsvertexbuffer.h"

#include <chemkit/constants.h>

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
    Q_UNUSED(type);

    buffer->draw();
}

void GraphicsPainter::drawSphere(GraphicsFloat radius)
{
    GraphicsSphere sphere(radius);
    GraphicsVertexBuffer *buffer = sphere.tesselate();

    draw(buffer);

    delete buffer;
}

void GraphicsPainter::drawSphere(const Point3g &center, GraphicsFloat radius)
{
    glPushMatrix();
    glTranslated(center.x(), center.y(), center.z());

    drawSphere(radius);

    glPopMatrix();
}

void GraphicsPainter::drawCylinder(GraphicsFloat radius, GraphicsFloat length)
{
    GraphicsCylinder cylinder(radius, length);
    GraphicsVertexBuffer *buffer = cylinder.tesselate(12, 10);

    draw(buffer);

    delete buffer;
}

void GraphicsPainter::drawCylinder(const Point3g &a, const Point3g &b, GraphicsFloat radius)
{
    glPushMatrix();

    glTranslatef(a.x(), a.y(), a.z());

    Vector3g vector = (a - b).normalized();
    Vector3g axis = vector.cross(-Vector3g::Z()).normalized();
    GraphicsFloat angle = vector.angle(-Vector3g::Z());
    glRotatef(-angle, axis.x(), axis.y(), axis.z());

    GraphicsFloat length = a.distance(b);

    drawCylinder(radius, length);

    glPopMatrix();
}

void GraphicsPainter::drawCircle(GraphicsFloat radius)
{
    Q_UNUSED(radius);
}

void GraphicsPainter::drawCircle(const Point3g &center, GraphicsFloat radius, const Vector3g &normal)
{
    Vector3g right(normal.y(), -normal.x(), 0); // vector orthogonal to normal

    glBegin(GL_TRIANGLE_FAN);

    glNormal3f(normal.x(), normal.y(), normal.z());
    glVertex3f(center.x(), center.y(), center.z());

    for(int angle = 0; angle <= 360; angle += 10){
        Point3g point = center.movedBy(radius, GraphicsQuaternion::rotate(right, normal, angle));

        glNormal3f(normal.x(), normal.y(), normal.z());
        glVertex3f(point.x(), point.y(), point.z());
    }

    glEnd();
}

void GraphicsPainter::drawTriangle(const Point3g &a, const Point3g &b, const Point3g &c)
{
    // verticies
    QVector<Point3g> verticies(3);
    verticies[0] = a;
    verticies[1] = b;
    verticies[2] = c;

    // normals
    Vector3g normal = (b - a).cross(c - b).normalized();
    QVector<Vector3g> normals(3);
    normals.fill(normal);

    // indicies
    QVector<unsigned short> indicies;
    indicies << 0 << 1 << 2;

    // setup buffer
    GraphicsVertexBuffer buffer;
    buffer.setVerticies(verticies);
    buffer.setNormals(normals);
    buffer.setIndicies(indicies);

    // draw buffer
    draw(&buffer);
}

void GraphicsPainter::drawRectangle(const Point3g &a, const Point3g &b, const Point3g &c, const Point3g &d)
{
    drawTriangle(a, b, c);
    drawTriangle(a, c, d);
}

void GraphicsPainter::drawSpline(const QList<Point3g> &points, GraphicsFloat radius, int order)
{
    if(points.size() < 2){
        return;
    }

    // calculate axis and up vectors (needed to calculate control points)
    QVector<Vector3g> axisVectors(points.size());
    axisVectors[0] = (points[1] - points[0]).normalized();

    QVector<Vector3g> upVectors(points.size());
    if(points.size() > 2)
        upVectors[0] = (points[1] - points[0]).cross(points[2] - points[1]).normalized();
    else
        upVectors[0] = Vector3g::Z();

    for(int i = 1; i < points.size(); i++){
        Vector3g axis = points[i] - points[i-1];

        if(i != (points.size() - 1)){
            GraphicsFloat angle = axis.angle(points[i+1] - points[i]);
            Vector3g rotationAxis = (points[i] - points[i-1]).cross(points[i+1] - points[i]).normalized();
            axis = GraphicsQuaternion::rotate(axis, rotationAxis, angle / 2.0);
        }

        axisVectors[i] = axis.normalized();

        Vector3g rotationAxis = axisVectors[i-1].cross(axisVectors[i]);
        GraphicsFloat angle = axis.angle(axisVectors[i-1]);
        Vector3g up = GraphicsQuaternion::rotate(upVectors[i-1], rotationAxis, angle);
        up.normalize();

        upVectors[i] = up;
    }

    // calculate control points
    QVector<Point3g> controlPoints(points.size() * 9);

    for(int i = 0; i < points.size(); i++){
        const Point3g &point = points.at(i);

        Vector3g upVector = upVectors[i];
        Vector3g rightVector = upVector.cross(axisVectors[i]);

        // 8 points around a square surrounding point

        // right
        Point3g right = point.movedBy(radius, rightVector);
        controlPoints[i*9+0] = right;

        // bottom right
        Point3g bottomRight = right.movedBy(-radius, upVector);
        controlPoints[i*9+1] = bottomRight;

        // bottom
        Point3g bottom = point.movedBy(-radius, upVector);
        controlPoints[i*9+2] = bottom;

        // bottom left
        Point3g bottomLeft = bottom.movedBy(-radius, rightVector);
        controlPoints[i*9+3] = bottomLeft;

        // left
        Point3g left = point.movedBy(-radius, rightVector);
        controlPoints[i*9+4] = left;

        // top left
        Point3g topLeft = left.movedBy(radius, upVector);
        controlPoints[i*9+5] = topLeft;

        // top
        Point3g top = point.movedBy(radius, upVector);
        controlPoints[i*9+6] = top;

        // top right
        Point3g topRight = top.movedBy(radius, rightVector);
        controlPoints[i*9+7] = topRight;

        // right (again)
        controlPoints[i*9+8] = right;
    }

    // build knot vector
    QVector<GraphicsFloat> uKnots(points.size() + order);
    for(int i = 0; i < uKnots.size(); i++){
        if(i < order)
            uKnots[i] = 0;
        else if(i > (points.size()) - 1)
            uKnots[i] = (points.size() + order) - 2 * order + 1;
        else
            uKnots[i] = (i - order + 1);
    }

    QVector<GraphicsFloat> vKnots(12);
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

void GraphicsPainter::drawNurbsSurface(const QVector<Point3g> &controlPoints, const QVector<GraphicsFloat> &uKnots, const QVector<GraphicsFloat> &vKnots, int uOrder, int vOrder)
{
    // nurbs control points
    QVector<GraphicsFloat> points(controlPoints.size() * 3);
    for(int i = 0; i < controlPoints.size(); i++){
        points[i*3+0] = controlPoints[i].x();
        points[i*3+1] = controlPoints[i].y();
        points[i*3+2] = controlPoints[i].z();
    }

    // nurbs knot vectors
    QVector<GraphicsFloat> uKnotVector(uKnots.size());
    for(int i = 0; i < uKnots.size(); i++)
        uKnotVector[i] = uKnots[i];

    QVector<GraphicsFloat> vKnotVector(vKnots.size());
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

    GraphicsFloat specular[] = {material->specularColor().redF(),
                                material->specularColor().greenF(),
                                material->specularColor().blueF(),
                                material->specularColor().alphaF()};
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
}

} // end chemkit namespace
