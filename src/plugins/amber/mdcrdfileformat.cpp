/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "mdcrdfileformat.h"

#include <boost/make_shared.hpp>

#include <chemkit/topology.h>
#include <chemkit/trajectory.h>
#include <chemkit/trajectoryfile.h>
#include <chemkit/trajectoryframe.h>
#include <chemkit/cartesiancoordinates.h>

MdcrdFileFormat::MdcrdFileFormat()
    : chemkit::TrajectoryFileFormat("mdcrd")
{
}

bool MdcrdFileFormat::read(std::istream &input, chemkit::TrajectoryFile *file)
{
    boost::shared_ptr<chemkit::Topology> topology = file->topology();
    if(!topology){
        setErrorString("Topology required to read 'mdcrd' trajectories.");
        return false;
    }

    boost::shared_ptr<chemkit::Trajectory> trajectory =
        boost::make_shared<chemkit::Trajectory>();

    trajectory->resize(topology->size());

    // comments line
    std::string comments;
    std::getline(input, comments);

    while(!input.eof()){
        chemkit::TrajectoryFrame *frame = trajectory->addFrame();

        for(size_t i = 0; i < trajectory->size(); i++){
            chemkit::Real x, y, z;
            input >> x >> y >> z;
            frame->setPosition(i, chemkit::Point3(x, y, z));
        }
    }

    file->setTrajectory(trajectory);

    return true;
}
