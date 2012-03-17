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

#include "trajectory.h"

#include <algorithm>

#include <chemkit/foreach.h>

#include "trajectoryframe.h"

namespace chemkit {

// === TrajectoryPrivate =================================================== //
class TrajectoryPrivate
{
public:
    size_t size;
    std::vector<TrajectoryFrame *> frames;
};

// === Trajectory ========================================================== //
/// \class Trajectory trajectory.h chemkit/trajectory.h
/// \ingroup chemkit-md
/// \brief The Trajectory class contains a molecular dynamics trajectory.
///
/// Trajectory objects contain multiple TrajectoryFrame objects. Each
/// trajectory frame contains the coordinates for each particle in the
/// system at a specific point in time.
///
/// Trajectories are usually associated with a Topology which contains
/// the atomic properties and atomic interactions for a system.
///
/// \see Topology, TrajectoryFrame, TrajectoryFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new trajectory with \p size.
Trajectory::Trajectory(size_t size)
    : d(new TrajectoryPrivate)
{
    d->size = size;
}

/// Destroys the trajectory object.
Trajectory::~Trajectory()
{
    foreach(TrajectoryFrame *frame, d->frames){
        delete frame;
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the number of particles in the trajectory to \p size.
void Trajectory::resize(size_t size)
{
    d->size = size;

    foreach(TrajectoryFrame *frame, d->frames){
        frame->resize(size);
    }
}

/// Returns the number of particles in the trajectory.
size_t Trajectory::size() const
{
    return d->size;
}

/// Returns \c true if the trajectory contains no frames.
bool Trajectory::isEmpty() const
{
    return frameCount() == 0;
}

// --- Frames -------------------------------------------------------------- //
/// Adds a new frame to the trajectory.
TrajectoryFrame* Trajectory::addFrame()
{
    TrajectoryFrame *frame = new TrajectoryFrame(this, d->size);
    d->frames.push_back(frame);
    return frame;
}

/// Removes \p frame from the trajectory.
bool Trajectory::removeFrame(TrajectoryFrame *frame)
{
    std::vector<TrajectoryFrame *>::iterator location = std::find(d->frames.begin(), d->frames.end(), frame);
    if(location == d->frames.end()){
        return false;
    }

    d->frames.erase(location);
    delete frame;

    return true;
}

/// Returns the frame at \p index in the trajectory.
TrajectoryFrame* Trajectory::frame(size_t index) const
{
    return d->frames[index];
}

/// Returns a list of the frames in the trajectory.
std::vector<TrajectoryFrame *> Trajectory::frames() const
{
    return d->frames;
}

/// Returns the number of frames in the trajectory.
size_t Trajectory::frameCount() const
{
    return d->frames.size();
}

} // end chemkit namespace
