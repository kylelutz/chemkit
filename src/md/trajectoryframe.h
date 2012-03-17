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

#ifndef CHEMKIT_TRAJECTORYFRAME_H
#define CHEMKIT_TRAJECTORYFRAME_H

#include "md.h"

#include <chemkit/point3.h>

namespace chemkit {

class UnitCell;
class CartesianCoordinates;
class Trajectory;
class TrajectoryFramePrivate;

class CHEMKIT_MD_EXPORT TrajectoryFrame
{
public:
    // properties
    size_t size() const;
    bool isEmpty() const;
    size_t index() const;
    Trajectory* trajectory() const;

    // time
    void setTime(Real time);
    Real time() const;

    // coordinates
    void setPosition(size_t index, const Point3 &position);
    Point3 position(size_t index) const;
    const CartesianCoordinates* coordinates() const;

    // unit cell
    void setUnitCell(UnitCell *cell);
    UnitCell* unitCell() const;

private:
    // construction and destruction
    TrajectoryFrame(Trajectory *trajectory, size_t size);
    ~TrajectoryFrame();

    void resize(size_t size);

    friend class Trajectory;

private:
    TrajectoryFramePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_TRAJECTORYFRAME_H
