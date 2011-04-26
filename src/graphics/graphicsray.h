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

#ifndef CHEMKIT_GRAPHICSRAY_H
#define CHEMKIT_GRAPHICSRAY_H

#include "graphics.h"

#include <chemkit/point3.h>
#include <chemkit/vector3.h>

namespace chemkit {

class CHEMKIT_GRAPHICS_EXPORT GraphicsRay
{
    public:
        // construction and destruction
        GraphicsRay();
        GraphicsRay(const Point3f &origin, const Vector3f &direction);
        GraphicsRay(const Point3f &origin, const Point3f &point);
        ~GraphicsRay();

        // properties
        void setOrigin(const Point3f &origin);
        Point3f origin() const;
        void setDirection(const Vector3f &direction);
        Vector3f direction() const;

        // geometry
        Point3f pointAt(float distance) const;
        bool intersectsSphere(const Point3f &center, float radius, float *distance = 0) const;
        bool intersectsCylinder(const Point3f &a, const Point3f &b, float radius, float *distance = 0) const;

    private:
        Point3f m_origin;
        Vector3f m_direction;
};

} // end chemkit namespace

#endif // CHEMKIT_GRAPHICSRAY_H
