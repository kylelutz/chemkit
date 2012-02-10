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

#ifndef CHEMKIT_GRAPHICSCAMERA_H
#define CHEMKIT_GRAPHICSCAMERA_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>

namespace chemkit {

class GraphicsView;
class GraphicsCameraPrivate;

class CHEMKIT_GRAPHICS_EXPORT GraphicsCamera
{
public:
    // construction and destruction
    GraphicsCamera();
    GraphicsCamera(const Point3f &position);
    GraphicsCamera(float x, float y, float z);
    ~GraphicsCamera();

    // position
    void setPosition(const Point3f &position);
    void setPosition(float x, float y, float z);
    Point3f position() const;
    float x() const;
    float y() const;
    float z() const;
    void moveTo(const Point3f &position);
    void moveTo(float x, float y, float z);
    void moveBy(const Vector3f &vector);
    void moveBy(float dx, float dy, float dz);
    void moveBy(float distance, const Vector3f &direction);
    void moveFoward(float distance);
    void moveBackward(float distance);
    void rotate(const Vector3f &axis, float angle, bool rotateDirection = true);
    void orbit(float dx, float dy, bool rotateDirection = true);
    void orbit(const Point3f &point, float dx, float dy, bool rotateDirection = true);

    // orientation
    void setDirection(const Vector3f &direction);
    Vector3f direction() const;
    void setFocus(const Point3f &point);
    Point3f focus() const;
    void lookAt(const Point3f &point);
    void setUpVector(const Vector3f &upVector);
    Vector3f upVector() const;
    void tilt(float angle);

private:
    // internal methods
    void setChanged(bool changed);
    bool changed() const;

    friend class GraphicsView;

private:
    GraphicsCameraPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSCAMERA_H
