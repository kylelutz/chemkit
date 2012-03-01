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

#include "grofileformat.h"

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/topology.h>
#include <chemkit/topologyfile.h>

GroFileFormat::GroFileFormat()
    : chemkit::TopologyFileFormat("gro")
{
}

GroFileFormat::~GroFileFormat()
{
}

bool GroFileFormat::read(std::istream &input, chemkit::TopologyFile *file)
{
    // read comments from the first line
    std::string comments;
    std::getline(input, comments);

    // read topology size from the second line
    size_t size = 0;
    std::string sizeString;
    std::getline(input, sizeString);
    boost::algorithm::trim(sizeString);

    try{
        size = boost::lexical_cast<size_t>(sizeString);
    }
    catch(boost::bad_lexical_cast&){
        setErrorString("Second line does not contain size.");
        return false;
    }

    boost::shared_ptr<chemkit::Topology> topology =
        boost::make_shared<chemkit::Topology>(size);

    for(size_t i = 0; i < size; i++){
        std::string line;
        std::getline(input, line);
        if(line.empty()){
            break;
        }

        boost::algorithm::trim_left(line);
        std::vector<std::string> tokens;
        boost::split(tokens,
                     line,
                     boost::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        if(tokens.size() < 2){
            break;
        }

        topology->setType(i, tokens[1]);
    }

    file->setTopology(topology);

    return true;
}

