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

#include "topologyfile.h"

#include <chemkit/topology.h>

namespace chemkit {

// === TopologyFilePrivate ================================================= //
class TopologyFilePrivate
{
public:
    boost::shared_ptr<Topology> topology;
};

// === TopologyFile ======================================================== //
/// \class TopologyFile topologyfile.h chemkit/topologyfile.h
/// \ingroup chemkit-md-io
/// \brief The TopologyFile class represents a file contatining a
///        molecular dynamics topology.
///
/// A list of supported topology file formats is available at:
/// http://wiki.chemkit.org/Features#Topology_File_Formats
///
/// \see Topology, TrajectoryFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new topology file.
TopologyFile::TopologyFile()
    : d(new TopologyFilePrivate)
{
}

/// Creates a new topology file with \p fileName.
TopologyFile::TopologyFile(const std::string &fileName)
    : GenericFile<TopologyFile, TopologyFileFormat>(fileName),
      d(new TopologyFilePrivate)
{
}

/// Destroys the topology file object.
TopologyFile::~TopologyFile()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns \c true if the file is empty.
bool TopologyFile::isEmpty() const
{
    return d->topology == 0;
}

// --- File Contents ------------------------------------------------------- //
/// Sets the topology to \p topology.
void TopologyFile::setTopology(const boost::shared_ptr<Topology> &topology)
{
    d->topology = topology;
}

/// Returns the topology in the file.
boost::shared_ptr<Topology> TopologyFile::topology() const
{
    return d->topology;
}

} // end chemkit namespace
