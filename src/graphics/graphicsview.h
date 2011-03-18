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

#ifndef CHEMKIT_GRAPHICSVIEW_H
#define CHEMKIT_GRAPHICSVIEW_H

#include "graphics.h"

#include "point3g.h"
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
        GraphicsView(GraphicsScene *scene, QWidget *parent = 0);
        ~GraphicsView();

        // properties
        void setScene(GraphicsScene *scene);
        GraphicsScene* scene() const;
        void setBackgroundColor(const QColor &color);
        QColor backgroundColor() const;
        void setTool(GraphicsTool *tool);
        GraphicsTool* tool() const;
        const GraphicsTransform& projectionTransform() const;
        const GraphicsTransform& modelViewTransform() const;

        // items
        void addItem(GraphicsItem *item);
        bool removeItem(GraphicsItem *item);
        bool deleteItem(GraphicsItem *item);
        QList<GraphicsItem *> items() const;
        int itemCount() const;

        // camera
        void setCamera(GraphicsCamera *camera);
        GraphicsCamera* camera() const;
        void setNearClipDistance(GraphicsFloat distance);
        GraphicsFloat nearClipDistance() const;
        void setFarClipDistance(GraphicsFloat distance);
        GraphicsFloat farClipDistance() const;
        QPointF project(const Point3g &point) const;
        Point3g unproject(qreal x, qreal y, qreal z) const;
        Point3g unproject(qreal x, qreal y, const Point3g &reference) const;
        GraphicsFloat depth(const Point3g &point) const;

        // lighting
        void addLight(GraphicsLight *light);
        bool removeLight(GraphicsLight *light);
        bool deleteLight(GraphicsLight *light);
        QList<GraphicsLight *> lights() const;
        int lightCount() const;
        GraphicsLight* light(int index = 0) const;

        // selection
        GraphicsItem* itemAt(int x, int y) const;
        QList<GraphicsItem *> itemsAt(int x, int y, bool sorted = true) const;

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
