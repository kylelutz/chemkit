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

#include "graphicsview.h"

#include <boost/make_shared.hpp>

#include <chemkit/foreach.h>
#include <chemkit/vector3.h>

#include "graphics.h"
#include "graphicsitem.h"
#include "graphicstool.h"
#include "graphicslight.h"
#include "graphicsscene.h"
#include "graphicscamera.h"
#include "graphicsoverlay.h"
#include "graphicspainter.h"
#include "graphicsmaterial.h"
#include "graphicstransform.h"
#include "graphicsnavigationtool.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE 0x809D
#endif

namespace chemkit {

// === GraphicsViewPrivate ================================================= //
class GraphicsViewPrivate
{
public:
    GraphicsViewPrivate();

    boost::shared_ptr<GraphicsScene> scene;
    boost::shared_ptr<GraphicsCamera> camera;
    boost::shared_ptr<GraphicsTool> tool;
    QColor backgroundColor;
    std::vector<boost::shared_ptr<GraphicsLight> > lights;
    GraphicsOverlay *overlay;
    bool overlayEnabled;
    GraphicsTransform modelViewTransform;
    GraphicsTransform projectionTransform;
    QGLShaderProgram *shader;
    float nearClipDistance;
    float farClipDistance;
    float fieldOfView;
    bool fogEnabled;
    bool hardwareIsSupported;
};

GraphicsViewPrivate::GraphicsViewPrivate()
{
    backgroundColor = Qt::black;
    overlay = new GraphicsOverlay;
    overlayEnabled = true;
    nearClipDistance = 0.01f;
    farClipDistance = 500.0f;
    fieldOfView = 45.0f;
    shader = 0;
    fogEnabled = true;
    hardwareIsSupported = false;
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
/// A gallery showing the different graphics items is available
/// at: http://wiki.chemkit.org/Graphics_Item_Gallery
///
/// The example below shows how to read a molecule from a file,
/// display it in a GraphicsView widget and setup mouse navigation.
/// \code
/// // read the molecule from the file
/// MoleculeFile file("/path/to/guanine.mol");
/// file.read();
/// boost::shared_ptr<Molecule> molecule = file.molecule();
///
/// // create the graphics view widget
/// GraphicsView view;
///
/// // add a molecule item to display the molecule
/// GraphicsMoleculeItem *item = new GraphicsMoleculeItem(molecule.get());
/// view.addItem(item);
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
    setScene(boost::make_shared<GraphicsScene>());
    setCamera(boost::make_shared<GraphicsCamera>(0, 0, 10));
    setTool(boost::make_shared<GraphicsNavigationTool>());
    setAutoFillBackground(false);
}

/// Creates a new graphics view widget displaying \p scene.
GraphicsView::GraphicsView(const boost::shared_ptr<GraphicsScene> &scene,
                           QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent),
      d(new GraphicsViewPrivate)
{
    setScene(scene);
    setCamera(boost::make_shared<GraphicsCamera>(0, 0, 10));
    setTool(boost::make_shared<GraphicsNavigationTool>());
    setAutoFillBackground(false);
}

/// Destroys the graphics view widget.
GraphicsView::~GraphicsView()
{
    if(d->scene){
        d->scene->removeView(this);
    }

    delete d->overlay;
    delete d->shader;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the graphics scene to show.
void GraphicsView::setScene(const boost::shared_ptr<GraphicsScene> &scene)
{
    if(scene == d->scene){
        return;
    }

    if(d->scene){
        d->scene->removeView(this);
    }

    d->scene = scene;

    if(d->scene){
        d->scene->addView(this);
    }

    update();
}

/// Returns the scene that the view is showing.
boost::shared_ptr<GraphicsScene> GraphicsView::scene() const
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
void GraphicsView::setTool(const boost::shared_ptr<GraphicsTool> &tool)
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
boost::shared_ptr<GraphicsTool> GraphicsView::tool() const
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
std::vector<GraphicsItem *> GraphicsView::items() const
{
    if(d->scene)
        return d->scene->items();

    return std::vector<GraphicsItem *>();
}

/// Returns the number of items in the view's scene.
///
/// \see GraphicsScene::itemCount()
size_t GraphicsView::itemCount() const
{
    if(d->scene)
        return d->scene->itemCount();

    return 0;
}

// --- Camera -------------------------------------------------------------- //
/// Sets the camera to \p camera.
void GraphicsView::setCamera(const boost::shared_ptr<GraphicsCamera> &camera)
{
    d->camera = camera;

    update();
}

/// Returns the camera.
boost::shared_ptr<GraphicsCamera> GraphicsView::camera() const
{
    return d->camera;
}

/// Projects a point from the scene to the window.
QPointF GraphicsView::project(const Point3f &point) const
{
    Eigen::Matrix<float, 4, 1> vector;
    vector[0] = point.x();
    vector[1] = point.y();
    vector[2] = point.z();
    vector[3] = 0;

    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    vector = transform.multiply(vector);
    vector *= 1.0 / vector[3];

    float winX = width() * (vector[0] + 1) / 2;
    float winY = height() * (vector[1] + 1) / 2;
    float winZ = (vector[2] + 1) / 2;

    // if winZ is greater than 1.0 the point is not
    // visible (it is either in front of the near clip
    // plane or behind the far clip plane).
    if(winZ > 1.0)
        return QPointF();
    else
        return QPointF(winX, height() - winY);
}

/// Unprojects a point from the window to the scene.
Point3f GraphicsView::unproject(qreal x, qreal y, qreal z) const
{
    // flip y
    y = height() - y;

    // adjust point to normalized window coordinates
    Eigen::Matrix<float, 4, 1> point;
    point[0] = 2 * x / width() - 1;
    point[1] = 2 * y / height() - 1;
    point[2] = 2 * z - 1;
    point[3] = 1;

    // map to object-space coordinates
    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    point = transform.inverseMultiply(point);
    point *= 1.0 / point[3];

    return Point3f(point[0], point[1], point[2]);
}

/// Unprojects a point from the window to the scene using the
/// depth of \p reference for the z coordinate.
Point3f GraphicsView::unproject(qreal x, qreal y, const Point3f &reference) const
{
    return unproject(x, y, depth(reference));
}

/// Returns the depth of point in the scene.
float GraphicsView::depth(const Point3f &point) const
{
    Eigen::Matrix<float, 4, 1> viewPoint;
    viewPoint[0] = point.x();
    viewPoint[1] = point.y();
    viewPoint[2] = point.z();
    viewPoint[3] = 1;

    GraphicsTransform transform = projectionTransform() * modelViewTransform();
    viewPoint = transform.multiply(viewPoint);
    viewPoint *= 1.0 / viewPoint[3];

    float winZ = (viewPoint[2] + 1) / 2;

    return winZ;
}

/// Sets the near clip distance.
void GraphicsView::setNearClipDistance(float distance)
{
    d->nearClipDistance = distance;
    update();
}

/// Returns the near clip distance.
float GraphicsView::nearClipDistance() const
{
    return d->nearClipDistance;
}

/// Sets the far clip distance.
void GraphicsView::setFarClipDistance(float distance)
{
    d->farClipDistance = distance;
    update();
}

/// Returns the far clip distance.
float GraphicsView::farClipDistance() const
{
    return d->farClipDistance;
}

// --- Lighting ------------------------------------------------------------ //
/// Adds \p light to the view.
void GraphicsView::addLight(const boost::shared_ptr<GraphicsLight> &light)
{
    d->lights.push_back(light);
}

/// Removes \p light from the view.
bool GraphicsView::removeLight(const boost::shared_ptr<GraphicsLight> &light)
{
    std::vector<boost::shared_ptr<GraphicsLight> >::iterator iter = std::find(d->lights.begin(),
                                                                              d->lights.end(),
                                                                              light);

    if(iter == d->lights.end()){
        return false;
    }

    d->lights.erase(iter);
    return true;
}

/// Returns a list of lights in the view.
std::vector<boost::shared_ptr<GraphicsLight> > GraphicsView::lights() const
{
    return d->lights;
}

/// Returns the number of lights in the view.
size_t GraphicsView::lightCount() const
{
    return d->lights.size();
}

/// Returns the light at \p index.
boost::shared_ptr<GraphicsLight> GraphicsView::light(size_t index) const
{
    assert(index < d->lights.size());

    return d->lights[index];
}

// --- Fog ----------------------------------------------------------------- //
/// Enables/disables for rendering for the view.
void GraphicsView::setFogEnabled(bool enabled)
{
    d->fogEnabled = enabled;
}

/// Returns \c true if fog is enabled.
bool GraphicsView::fogEnabled() const
{
    return d->fogEnabled;
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
std::vector<GraphicsItem *> GraphicsView::itemsAt(int x, int y, bool sorted) const
{
    if(!d->scene){
        return std::vector<GraphicsItem *>();
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
        qWarning() << "Error: OpenGL version is not 2.0 or later.";
        d->hardwareIsSupported = false;
        return;
    }

    d->hardwareIsSupported = true;

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
    if(!d->scene || !d->camera){
        return;
    }

    d->shader->bind();

    // clear
    qglClearColor(d->backgroundColor);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // setup camera (if necessary)
    if(camera()->changed()){
        glMatrixMode(GL_MODELVIEW);

        Vector3f f = camera()->direction();
        Vector3f s = f.cross(camera()->upVector());
        Vector3f u = s.cross(f);

        Eigen::Matrix<float, 4, 4> transform;
        transform <<  s.x(),  s.y(),  s.z(), 0.0f,
                      u.x(),  u.y(),  u.z(), 0.0f,
                     -f.x(), -f.y(), -f.z(), 0.0f,
                       0.0f,   0.0f,   0.0f, 1.0f;
        d->modelViewTransform = GraphicsTransform(transform);

        d->modelViewTransform *= GraphicsTransform::translation(-camera()->position());
        glLoadMatrixf(d->modelViewTransform.data());

        camera()->setChanged(false);
    }

    // setup fog
    float fogColor[] = {d->backgroundColor.redF(),
                        d->backgroundColor.greenF(),
                        d->backgroundColor.blueF()};
    glFogfv(GL_FOG_COLOR, fogColor);

    if(d->fogEnabled){
        glFogf(GL_FOG_START, nearClipDistance());
        glFogf(GL_FOG_END, farClipDistance());
    }
    else{
        glFogf(GL_FOG_START, farClipDistance());
        glFogf(GL_FOG_END, farClipDistance());
    }

    // draw items
    GraphicsPainter painter;

    std::vector<GraphicsItem *> nonOpaqueItems;

    foreach(GraphicsItem *item, scene()->items()){
        if(!item->isVisible())
            continue;

        if(!item->isOpaque()){
            nonOpaqueItems.push_back(item);
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

    if(!nonOpaqueItems.empty()){
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
    float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    d->projectionTransform = GraphicsTransform::perspective(d->fieldOfView,
                                                            aspectRatio,
                                                            d->nearClipDistance,
                                                            d->farClipDistance);
    glLoadMatrixf(d->projectionTransform.data());
    glMatrixMode(GL_MODELVIEW);

    // resize overlay
    d->overlay->setSceneRect(0, 0, width, height);

    camera()->setChanged(true);
}

// --- Events -------------------------------------------------------------- //
void GraphicsView::paintEvent(QPaintEvent *event)
{
    if(!d->hardwareIsSupported){
        // if the graphics hardware doesn't support OpenGL 2.0 or
        // later just draw an error message on screen and return
        QPainter painter(this);
        painter.setPen(Qt::white);
        painter.setBrush(Qt::white);
        painter.drawText(QPointF(5, 25),
                         "Error: OpenGL 2.0 not supported by hardware.");
        return;
    }

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
    Point3f nearPoint = unproject(x, y, 0);
    Point3f farPoint = unproject(x, y, 1);

    return GraphicsRay(nearPoint, farPoint);
}

} // end chemkit namespace
