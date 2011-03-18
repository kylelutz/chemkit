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

#include <QtTest>

#include <chemkit/graphicsview.h>
#include <chemkit/graphicscamera.h>

class GraphicsCameraTest : public QObject
{
    Q_OBJECT

    private slots:
        void view();
        void position();
        void direction();
        void upVector();
};

void GraphicsCameraTest::view()
{
    chemkit::GraphicsCamera camera;
    QVERIFY(camera.view() == 0);

    chemkit::GraphicsView view;
    view.setCamera(&camera);
    QVERIFY(camera.view() == &view);

    view.setCamera(0);
    QVERIFY(camera.view() == 0);
}

void GraphicsCameraTest::position()
{
    chemkit::GraphicsCamera camera;
    QCOMPARE(camera.position(), chemkit::Point3g(0, 0, 0));

    camera.setPosition(10, 11, 12);
    QCOMPARE(camera.position(), chemkit::Point3g(10, 11, 12));
    QCOMPARE(camera.x(), chemkit::GraphicsFloat(10));
    QCOMPARE(camera.y(), chemkit::GraphicsFloat(11));
    QCOMPARE(camera.z(), chemkit::GraphicsFloat(12));

    camera.setPosition(chemkit::Point3g(1, -2, 3));
    QCOMPARE(camera.position(), chemkit::Point3g(1, -2, 3));
}

void GraphicsCameraTest::direction()
{
    chemkit::GraphicsCamera camera;

    // default direction is on the -Z axis
    QCOMPARE(camera.direction(), chemkit::Vector3g(0, 0, -1));

    camera.setDirection(chemkit::Vector3g(1, 0, 0));
    QCOMPARE(camera.direction(), chemkit::Vector3g(1, 0, 0));

    // ensure direction is normalized when set
    camera.setDirection(chemkit::Vector3g(0, 5, 0));
    QCOMPARE(camera.direction(), chemkit::Vector3g(0, 1, 0));
}

void GraphicsCameraTest::upVector()
{
    chemkit::GraphicsCamera camera;

    // default up vector is on the Y axis
    QCOMPARE(camera.upVector(), chemkit::Vector3g(0, 1, 0));

    camera.setUpVector(chemkit::Vector3g(1, 0, 0));
    QCOMPARE(camera.upVector(), chemkit::Vector3g(1, 0, 0));

    // ensure up vector is normalized when set
    camera.setUpVector(chemkit::Vector3g(-3, 0, 0));
    QCOMPARE(camera.upVector(), chemkit::Vector3g(-1, 0, 0));
}

QTEST_MAIN(GraphicsCameraTest)
#include "graphicscameratest.moc"
