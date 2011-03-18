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

#include "graphicsview.h"

#include "graphics.h"
#include "graphicsitem.h"
#include "graphicstool.h"
#include "graphicslight.h"
#include "graphicsscene.h"
#include "graphicscamera.h"
#include "graphicsvector.h"
#include "graphicsoverlay.h"
#include "graphicspainter.h"
#include "graphicsmaterial.h"
#include "graphicstransform.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE 0x809D
#endif

namespace chemkit {

// === GraphicsViewPrivate ================================================= //
class GraphicsViewPrivate
{
    public:
        GraphicsViewPrivate();

        GraphicsScene *scene;
        bool ownScene;
        GraphicsCamera *camera;
        bool ownCamera;
        GraphicsTool *tool;
        QColor backgroundColor;
        QList<GraphicsLight *> lights;
        GraphicsOverlay *overlay;
        bool overlayEnabled;
        GraphicsTransform modelViewTransform;
        GraphicsTransform projectionTransform;
        QGLShaderProgram *shader;
        GraphicsFloat nearClipDistance;
        GraphicsFloat farClipDistance;
        GraphicsFloat fieldOfView;
};

GraphicsViewPrivate::GraphicsViewPrivate()
{
    scene = 0;
    tool = 0;
    backgroundColor = Qt::black;
    camera = 0;
    overlay = new GraphicsOverlay;
    overlayEnabled = true;
    nearClipDistance = 0.01f;
    farClipDistance = 500.0f;
    fieldOfView = 45.0f;
    shader = 0;
}

// === GraphicsView ======================================================== //
/// \class GraphicsView graphicsview.h chemkit/graphicsview.h
/// \ingroup chemkit-graphics
/// \brief The GraphicsView class provides a widget for molecular
///        visualization.
///
/// The GraphicsView class is the central component of the
/// chemkit-graphics library and is responsible for displaying all
/// types of graphical objects.
///
/// The following classes are used to display molecules, proteins,
/// and nucleic acids:
///   - GraphicsMoleculeItem
///   - GraphicsProteinItem
///   - GraphicsNucleicAcidItem
///
/// The example below shows how to read a molecule from a file,
/// display it in a GraphicsView widget and setup mouse navigation.
/// \code
/// // read the molecule from the file
/// ChemicalFile file("/path/to/guanine.mol");
/// file.read();
/// Molecule *molecule = file.molecule();
///
/// // create the graphics view widget
/// GraphicsView view;
///
/// // add a molecule item to display the molecule
/// GraphicsMoleculeItem *item = new GraphicsMoleculeItem(molecule);
/// view.addItem(item);
///
/// // add the navigation tool for mouse navigation
/// view.setTool(new GraphicsNavigationTool);
///
/// // run the application
/// QApplication app;
/// view.show();
/// app.exec();
/// \endcode
///
/// The result will be a window showing the guanine molecule:
/// \image html graphics-view-example.png

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new graphics view widget.
GraphicsView::GraphicsView(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent),
      d(new GraphicsViewPrivate)
{
    setScene(new GraphicsScene);
    d->ownScene = true;
    setCamera(new GraphicsCamera(0, 0, 10));
    d->ownCamera = true;
    setAutoFillBackground(false);
}

/// Creates a new graphics view widget displaying \p scene.
GraphicsView::GraphicsView(GraphicsScene *scene, QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent),
      d(new GraphicsViewPrivate)
{
    setScene(scene);
    d->ownScene = false;
    setCamera(new GraphicsCamera(0, 0, 10));
    d->ownCamera = true;
    setAutoFillBackground(false);
}

/// Destroys the graphics view widget.
GraphicsView::~GraphicsView()
{
    if(d->scene){
        d->scene->removeView(this);

        if(d->ownScene){
            delete d->scene;
        }
    }

    if(d->camera){
        d->camera->setView(0);

        if(d->ownCamera){
            delete d->camera;
        }
    }

    delete d->tool;
    delete d->overlay;
    delete d->shader;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the graphics scene to show.
void GraphicsView::setScene(GraphicsScene *scene)
{
    if(scene == d->scene){
        return;
    }

    if(d->scene){
        d->scene->removeView(this);

        if(d->ownScene){
            delete d->scene;
        }
    }

    d->scene = scene;
    d->ownScene = false;

    if(d->scene){
        d->scene->addView(this);
    }

    update();
}

/// Returns the scene that the view is showing.
GraphicsScene* GraphicsView::scene() const
{
    return d->scene;
}

/// Sets the background color.
void GraphicsView::setBackgroundColor(const QColor &color)
{
    d->backgroundColor = color;
    update();
}

/// Returns the background color.
QColor GraphicsView::backgroundColor() const
{
    return d->backgroundColor;
}

/// Sets the current tool to \p tool.
///
/// The view takes ownership of the tool.
void GraphicsView::setTool(GraphicsTool *tool)
{
    if(d->tool){
        // notify the current tool that the tool changed
        d->tool->toolChanged(tool);
    }

    d->tool = tool;

    if(d->tool){
        d->tool->setView(this);

        // notify the new tool that the tool changed
        d->tool->toolChanged(tool);
    }
}

/// Returns the current tool.
GraphicsTool* GraphicsView::tool() const
{
    return d->tool;
}

/// Returns the projection transformation.
const GraphicsTransform& GraphicsView::projectionTransform() const
{
    return d->projectionTransform;
}

/// Returns the model view transformation.
const GraphicsTransform& GraphicsView::modelViewTransform() const
{
    return d->modelViewTransform;
}

// --- Items --------------------------------------------------------------- //
/// Adds \p item to the view's scene.
///
/// The scene takes ownership of the item.
///
/// Equivalent to:
/// \code
/// scene()->addItem(item);
/// \endcode
///
/// \see GraphicsScene::addItem()
void GraphicsView::addItem(GraphicsItem *item)
{
    if(d->scene)
        d->scene->addItem(item);
}

/// Removes \p item from the view's scene. Returns \c true if the
/// item was found and removed successfully.
///
/// The ownership of the item is passed to the caller.
///
/// Equivalent to:
/// \code
/// scene()->removeItem(item);
/// \endcode
///
/// \see GraphicsScene::removeItem()
bool GraphicsView::removeItem(GraphicsItem *item)
{
    if(d->scene)
        return d->scene->removeItem(item);

    return false;
}

/// Removes \p item from the view's scene and deletes it. Returns
/// \c true if the item was found and deleted successfully.
///
/// Equivalent to:
/// \code
/// scene()->deleteItem(item);
/// \endcode
///
/// \see GraphicsScene::deleteItem()
bool GraphicsView::deleteItem(GraphicsItem *item)
{
    if(d->scene)
        return d->scene->deleteItem(item);

    return false;
}

/// Returns a list of all the items in the view's scene.
///
/// \see GraphicsScene::items()
QList<GraphicsItem *> GraphicsView::items() const
{
    if(d->scene)
        return d->scene->items();

    return QList<GraphicsItem *>();
}

/// Returns the number of items in the view's scene.
///
/// \see GraphicsScene::itemCount()
int GraphicsView::itemCount() const
{
    if(d->scene)
        return d->scene->itemCount();

    return 0;
}

// --- Camera -------------------------------------------------------------- //
/// Sets the camera to \p camera.
///
/// The view takes ownership of the camera.
void GraphicsView::setCamera(GraphicsCamera *camera)
{
    if(d->camera){
        d->camera->setView(0);

        if(d->ownCamera){
            delete d->camera;
        }
    }

    d->camera = camera;
    d->ownCamera = false;

    if(d->camera){
        camera->setView(this);
    }

    update();
}

/// Returns the camera.
GraphicsCamera* GraphicsView::camera() const
{
    return d->camera;
}

/// Projects a point from the scene to the window.
QPointF GraphicsView::project(const Point3g &point) const
{
    StaticVector<GraphicsFloat, 4> vector(point.data(), 3);
    vector[3] = 0;

    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    vector = transform.multiply(vector);
    vector.scale(1.0 / vector[3]);

    GraphicsFloat winX = width() * (vector[0] + 1) / 2;
    GraphicsFloat winY = height() * (vector[1] + 1) / 2;
    GraphicsFloat winZ = (vector[2] + 1) / 2;

    // if winZ is greater than 1.0 the point is not
    // visible (it is either in front of the near clip
    // plane or behind the far clip plane).
    if(winZ > 1.0)
        return QPointF();
    else
        return QPointF(winX, height() - winY);
}

/// Unprojects a point from the window to the scene.
Point3g GraphicsView::unproject(qreal x, qreal y, qreal z) const
{
    // flip y
    y = height() - y;

    // adjust point to normalized window coordinates
    StaticVector<GraphicsFloat, 4> point;
    point[0] = 2 * x / width() - 1;
    point[1] = 2 * y / height() - 1;
    point[2] = 2 * z - 1;
    point[3] = 1;

    // map to object-space coordinates
    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    point = transform.inverseMultiply(point);
    point.scale(1.0 / point.value(3));

    return Point3g(point[0], point[1], point[2]);
}

/// Unprojects a point from the window to the scene using the
/// depth of \p reference for the z coordinate.
Point3g GraphicsView::unproject(qreal x, qreal y, const Point3g &reference) const
{
    return unproject(x, y, depth(reference));
}

/// Returns the depth of point in the scene.
GraphicsFloat GraphicsView::depth(const Point3g &point) const
{
    StaticVector<GraphicsFloat, 4> viewPoint(point.data(), 3);
    viewPoint[3] = 1;

    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    viewPoint = transform.multiply(viewPoint);
    viewPoint.scale(1.0 / viewPoint[3]);

    GraphicsFloat winZ = (viewPoint[2] + 1) / 2;

    return winZ;
}

/// Sets the near clip distance.
void GraphicsView::setNearClipDistance(GraphicsFloat distance)
{
    d->nearClipDistance = distance;
    update();
}

/// Returns the near clip distance.
GraphicsFloat GraphicsView::nearClipDistance() const
{
    return d->nearClipDistance;
}

/// Sets the far clip distance.
void GraphicsView::setFarClipDistance(GraphicsFloat distance)
{
    d->farClipDistance = distance;
    update();
}

/// Returns the far clip distance.
GraphicsFloat GraphicsView::farClipDistance() const
{
    return d->farClipDistance;
}

// --- Lighting ------------------------------------------------------------ //
/// Adds \p light to the view.
///
/// The view takes ownership of the light.
void GraphicsView::addLight(GraphicsLight *light)
{
    d->lights.append(light);
}

/// Removes \p light from the view.
///
/// The ownership of the light is passed to the caller.
bool GraphicsView::removeLight(GraphicsLight *light)
{
    return d->lights.removeOne(light);
}

/// Removes \p light from the view and deletes it.
///
/// Equivalent to:
/// \code
/// view.removeLight(light);
/// delete light;
/// \endcode
bool GraphicsView::deleteLight(GraphicsLight *light)
{
    bool found = removeLight(light);

    if(found){
        delete light;
    }

    return found;
}

/// Returns a list of lights in the view.
QList<GraphicsLight *> GraphicsView::lights() const
{
    return d->lights;
}

/// Returns the number of lights in the view.
int GraphicsView::lightCount() const
{
    return lights().size();
}

/// Returns the light at \p index.
GraphicsLight* GraphicsView::light(int index) const
{
    return d->lights.value(index, 0);
}

// --- Selection ----------------------------------------------------------- //
/// Returns the item at the window position (\p x, \p y).
GraphicsItem* GraphicsView::itemAt(int x, int y) const
{
    if(!d->scene){
        return 0;
    }

    GraphicsRay ray = buildPickRay(x, y);

    return scene()->item(ray);
}

/// Returns a list of all items under the window point (\p x, \p y).
QList<GraphicsItem *> GraphicsView::itemsAt(int x, int y, bool sorted) const
{
    if(!d->scene){
        return QList<GraphicsItem *>();
    }

    GraphicsRay ray = buildPickRay(x, y);

    return scene()->items(ray, sorted);
}

// --- Overlay ------------------------------------------------------------- //
/// Returns the overlay for the scene.
GraphicsOverlay* GraphicsView::overlay() const
{
    return d->overlay;
}

/// Sets whether or not the overlay is enabled.
void GraphicsView::setOverlayEnabled(bool enabled)
{
    d->overlayEnabled = enabled;
}

/// Returns \c true if the overlay is enabled.
bool GraphicsView::overlayEnabled() const
{
    return d->overlayEnabled;
}

// --- OpenGL -------------------------------------------------------------- //
void GraphicsView::initializeGL()
{
    // check opengl version
    if(!(QGLFormat::openGLVersionFlags() & QGLFormat::OpenGL_Version_2_0)){
        qWarning() << "GraphicsView: OpenGL version is not 2.0 or later.";
    }

    // background color
    qglClearColor(d->backgroundColor);

    // opengl settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_MULTISAMPLE);

    // materials
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // alpha blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    d->shader = new QGLShaderProgram(this);
    bool ok;
    ok = d->shader->addShaderFromSourceFile(QGLShader::Vertex, ":/shaders/phong.vert");
    if(!ok){
        qDebug() << "GraphicsView: Failed to load vertex shader:" << d->shader->log();
    }

    ok = d->shader->addShaderFromSourceFile(QGLShader::Fragment, ":/shaders/phong.frag");
    if(!ok){
        qDebug() << "GraphicsView: Failed to load fragment shader" << d->shader->log();
    }
}

void GraphicsView::paintGL()
{
    if(!d->scene){
        return;
    }

    d->shader->bind();

    // clear
    qglClearColor(d->backgroundColor);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // setup camera (if necessary)
    if(camera()->changed()){
        glMatrixMode(GL_MODELVIEW);

        GraphicsVector f = camera()->direction();
        GraphicsVector s = f.cross(camera()->upVector());
        GraphicsVector u = s.cross(f);

        d->modelViewTransform << s.x(),  s.y(),  s.z(), 0,
                                 u.x(),  u.y(),  u.z(), 0,
                                -f.x(), -f.y(), -f.z(), 0,
                                     0,      0,      0, 1;

        d->modelViewTransform *= GraphicsTransform::translation(-camera()->position());
        glLoadMatrixf(d->modelViewTransform.data());

        camera()->setChanged(false);
    }

    // draw items
    GraphicsPainter painter;

    QList<GraphicsItem *> nonOpaqueItems;

    foreach(GraphicsItem *item, scene()->items()){
        if(!item->isVisible())
            continue;

        if(!item->isOpaque()){
            nonOpaqueItems.append(item);
        }
        else{
            glPushMatrix();

            GraphicsTransform transform = item->transform();
            glMultMatrixf(transform.data());

            painter.setMaterial(item->material());
            item->paint(&painter);

            glPopMatrix();
        }
    }

    if(!nonOpaqueItems.isEmpty()){
        glEnable(GL_BLEND);

        foreach(GraphicsItem *item, nonOpaqueItems){
            glPushMatrix();

            GraphicsTransform transform = item->transform();
            glMultMatrixf(transform.data());

            painter.setMaterial(item->material());
            item->paint(&painter);

            glPopMatrix();
        }

        glDisable(GL_BLEND);
    }
}

void GraphicsView::resizeGL(int width, int height)
{
    // resize viewport
    glViewport(0, 0, width, height);

    // setup projection matrix
    glMatrixMode(GL_PROJECTION);
    GraphicsFloat aspectRatio = static_cast<GraphicsFloat>(width) / static_cast<GraphicsFloat>(height);
    d->projectionTransform = GraphicsTransform::perspective(d->fieldOfView, aspectRatio, d->nearClipDistance, d->farClipDistance);
    glLoadMatrixf(d->projectionTransform.data());
    glMatrixMode(GL_MODELVIEW);

    // resize overlay
    d->overlay->setSceneRect(0, 0, width, height);

    camera()->setChanged(true);
}

// --- Events -------------------------------------------------------------- //
void GraphicsView::paintEvent(QPaintEvent *event)
{
    // draw opengl
    makeCurrent();
    paintGL();

    // draw overlay
    if(overlayEnabled()){
        d->overlay->updateBindings(this);

        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);

        QPainter painter(this);
        painter.setBackgroundMode(Qt::TransparentMode);
        painter.setRenderHint(QPainter::Antialiasing);
        painter.setRenderHint(QPainter::TextAntialiasing);
        d->overlay->render(&painter);
        painter.end();

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
    }

    event->accept();
}

void GraphicsView::mousePressEvent(QMouseEvent *event)
{
    if(d->tool)
        d->tool->mousePressEvent(event);
}

void GraphicsView::mouseReleaseEvent(QMouseEvent *event)
{
    if(d->tool)
        d->tool->mouseReleaseEvent(event);
}

void GraphicsView::mouseDoubleClickEvent(QMouseEvent *event)
{
    if(d->tool)
        d->tool->mouseDoubleClickEvent(event);
}

void GraphicsView::mouseMoveEvent(QMouseEvent *event)
{
    if(d->tool)
        d->tool->mouseMoveEvent(event);
}

void GraphicsView::wheelEvent(QWheelEvent *event)
{
    if(d->tool)
        d->tool->wheelEvent(event);
}

// --- Internal Methods ---------------------------------------------------- //
GraphicsRay GraphicsView::buildPickRay(int x, int y) const
{
    Point3g nearPoint = unproject(x, y, 0);
    Point3g farPoint = unproject(x, y, 1);

    return GraphicsRay(nearPoint, farPoint);
}

} // end chemkit namespace
