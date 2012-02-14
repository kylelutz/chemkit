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

#ifndef CHEMKIT_GRAPHICSVIEW_H
#define CHEMKIT_GRAPHICSVIEW_H

#include "graphics.h"

#include <boost/shared_ptr.hpp>

#include <chemkit/point3.h>

#include "graphicsray.h"

namespace chemkit {

class GraphicsItem;
class GraphicsTool;
class GraphicsLight;
class GraphicsScene;
class GraphicsCamera;
class GraphicsOverlay;
class GraphicsTransform;
class GraphicsViewPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsView : public QGLWidget
{
    Q_OBJECT

public:
    // construction and destruction
    GraphicsView(QWidget *parent = 0);
    GraphicsView(const boost::shared_ptr<GraphicsScene> &scene, QWidget *parent = 0);
    ~GraphicsView();

    // properties
    void setScene(const boost::shared_ptr<GraphicsScene> &scene);
    boost::shared_ptr<GraphicsScene> scene() const;
    void setBackgroundColor(const QColor &color);
    QColor backgroundColor() const;
    void setTool(const boost::shared_ptr<GraphicsTool> &tool);
    boost::shared_ptr<GraphicsTool> tool() const;
    const GraphicsTransform& projectionTransform() const;
    const GraphicsTransform& modelViewTransform() const;

    // items
    void addItem(GraphicsItem *item);
    bool removeItem(GraphicsItem *item);
    bool deleteItem(GraphicsItem *item);
    std::vector<GraphicsItem *> items() const;
    size_t itemCount() const;

    // camera
    void setCamera(const boost::shared_ptr<GraphicsCamera> &camera);
    boost::shared_ptr<GraphicsCamera> camera() const;
    void setNearClipDistance(float distance);
    float nearClipDistance() const;
    void setFarClipDistance(float distance);
    float farClipDistance() const;
    QPointF project(const Point3f &point) const;
    Point3f unproject(qreal x, qreal y, qreal z) const;
    Point3f unproject(qreal x, qreal y, const Point3f &reference) const;
    float depth(const Point3f &point) const;

    // lighting
    void addLight(const boost::shared_ptr<GraphicsLight> &light);
    bool removeLight(const boost::shared_ptr<GraphicsLight> &light);
    std::vector<boost::shared_ptr<GraphicsLight> > lights() const;
    size_t lightCount() const;
    boost::shared_ptr<GraphicsLight> light(size_t index = 0) const;

    // fog
    void setFogEnabled(bool enabled);
    bool fogEnabled() const;

    // selection
    GraphicsItem* itemAt(int x, int y) const;
    std::vector<GraphicsItem *> itemsAt(int x, int y, bool sorted = true) const;

    // overlay
    GraphicsOverlay* overlay() const;
    void setOverlayEnabled(bool enabled);
    bool overlayEnabled() const;

protected:
    // opengl
    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int width, int height);

    // events
    virtual void paintEvent(QPaintEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseDoubleClickEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent *event);

private:
    // internal methods
    GraphicsRay buildPickRay(int x, int y) const;

private:
    GraphicsViewPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSVIEW_H
