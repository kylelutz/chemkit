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

#include "graphicscameratest.h"

#include <chemkit/graphicsview.h>
#include <chemkit/graphicscamera.h>

void GraphicsCameraTest::position()
{
    chemkit::GraphicsCamera camera;
    QCOMPARE(camera.position(), chemkit::Point3f(0, 0, 0));

    camera.setPosition(10, 11, 12);
    QCOMPARE(camera.position(), chemkit::Point3f(10, 11, 12));
    QCOMPARE(camera.x(), 10.0f);
    QCOMPARE(camera.y(), 11.0f);
    QCOMPARE(camera.z(), 12.0f);

    camera.setPosition(chemkit::Point3f(1, -2, 3));
    QCOMPARE(camera.position(), chemkit::Point3f(1, -2, 3));
}

void GraphicsCameraTest::direction()
{
    chemkit::GraphicsCamera camera;

    // default direction is on the -Z axis
    QCOMPARE(camera.direction(), chemkit::Vector3f(0, 0, -1));

    camera.setDirection(chemkit::Vector3f(1, 0, 0));
    QCOMPARE(camera.direction(), chemkit::Vector3f(1, 0, 0));

    // ensure direction is normalized when set
    camera.setDirection(chemkit::Vector3f(0, 5, 0));
    QCOMPARE(camera.direction(), chemkit::Vector3f(0, 1, 0));
}

void GraphicsCameraTest::upVector()
{
    chemkit::GraphicsCamera camera;

    // default up vector is on the Y axis
    QCOMPARE(camera.upVector(), chemkit::Vector3f(0, 1, 0));

    camera.setUpVector(chemkit::Vector3f(1, 0, 0));
    QCOMPARE(camera.upVector(), chemkit::Vector3f(1, 0, 0));

    // ensure up vector is normalized when set
    camera.setUpVector(chemkit::Vector3f(-3, 0, 0));
    QCOMPARE(camera.upVector(), chemkit::Vector3f(-1, 0, 0));
}

QTEST_MAIN(GraphicsCameraTest)
