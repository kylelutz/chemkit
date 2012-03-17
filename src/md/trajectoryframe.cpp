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

#include "trajectoryframe.h"

#include <vector>
#include <algorithm>

#include <chemkit/unitcell.h>
#include <chemkit/cartesiancoordinates.h>

#include "trajectory.h"

namespace chemkit {

// === TrajectoryFramePrivate ============================================== //
class TrajectoryFramePrivate
{
public:
    Trajectory *trajectory;
    Real time;
    CartesianCoordinates *coordinates;
    UnitCell *unitCell;
};

// === TrajectoryFrame ===================================================== //
/// \class TrajectoryFrame trajectoryframe.h chemkit/trajectoryframe.h
/// \ingroup chemkit-md
/// \brief The TrajectoryFrame class represents a single frame in a
///        trajectory.
///
/// TrajectoryFrame objects are created with the
/// Trajectory::addFrame() method and destroyed with the
/// Trajectory::removeFrame() method.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new trajectory frame.
TrajectoryFrame::TrajectoryFrame(Trajectory *trajectory, size_t size)
    : d(new TrajectoryFramePrivate)
{
    d->trajectory = trajectory;
    d->time = 0;
    d->coordinates = new CartesianCoordinates(size);
    d->unitCell = 0;
}

/// Destroys the trajectory frame object.
TrajectoryFrame::~TrajectoryFrame()
{
    delete d->coordinates;
    delete d->unitCell;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the number of coordinates in the frame to \p size.
void TrajectoryFrame::resize(size_t size)
{
    d->coordinates->resize(size);
}

/// Returns the number of coordinates in the frame.
size_t TrajectoryFrame::size() const
{
    return d->coordinates->size();
}

/// Returns \c true if the frame contains no coordinates.
bool TrajectoryFrame::isEmpty() const
{
    return d->coordinates->isEmpty();
}

/// Returns the index of the frame in the trajectory.
size_t TrajectoryFrame::index() const
{
    const std::vector<TrajectoryFrame *> &frames = trajectory()->frames();

    return std::distance(frames.begin(), std::find(frames.begin(), frames.end(), this));
}

/// Returns the trajectory that the frame belongs to.
Trajectory* TrajectoryFrame::trajectory() const
{
    return d->trajectory;
}

// --- Time ---------------------------------------------------------------- //
/// Sets the time for the trajectory frame to \p time.
void TrajectoryFrame::setTime(Real time)
{
    d->time = time;
}

/// Returns the time of the trajectory frame.
Real TrajectoryFrame::time() const
{
    return d->time;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the coordinates at \p index to \p position.
void TrajectoryFrame::setPosition(size_t index, const Point3 &position)
{
    d->coordinates->setPosition(index, position);
}

/// Returns the position at \p index.
Point3 TrajectoryFrame::position(size_t index) const
{
    return d->coordinates->position(index);
}

/// Returns the coordinates for the frame.
const CartesianCoordinates* TrajectoryFrame::coordinates() const
{
    return d->coordinates;
}

// --- Unit Cell ----------------------------------------------------------- //
/// Sets the unit cell for the frame to \p cell.
void TrajectoryFrame::setUnitCell(UnitCell *cell)
{
    d->unitCell = cell;
}

/// Returns the unit cell for the frame.
UnitCell* TrajectoryFrame::unitCell() const
{
    return d->unitCell;
}

} // end chemkit namespace
