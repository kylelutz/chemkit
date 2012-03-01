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
** OWNER OR CONTRIBUTORS BE LIABLgE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include "topfileformat.h"

#include <boost/make_shared.hpp>

#include <chemkit/topology.h>
#include <chemkit/topologyfile.h>

TopFileFormat::TopFileFormat()
    : chemkit::TopologyFileFormat("top")
{
}

TopFileFormat::~TopFileFormat()
{
}

bool TopFileFormat::read(std::istream &input, chemkit::TopologyFile *file)
{
    boost::shared_ptr<chemkit::Topology> topology =
        boost::make_shared<chemkit::Topology>();

    enum Section {
        Unknown,
        MoleculeType,
        Atoms,
        Bonds,
        Pairs,
        Angles,
        Dihedrals,
        PositionRestraints,
        System,
        Molecules
    };

    Section section = Unknown;

    while(!input.eof()){
        std::string line;
        std::getline(input, line);

        if(line.empty()){
            continue; // blank line
        }
        else if(line[0] == ';'){
            continue; // comment line
        }
        else if(line[0] == '#'){
            continue; // preprocessor command
        }
        else if(line[0] == '['){
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            std::string sectionName = tokens[1];

            if(sectionName == "moleculetype"){
                section = MoleculeType;
            }
            else if(sectionName == "atoms"){
                section = Atoms;
            }
            else if(sectionName == "bonds"){
                section = Bonds;
            }
            else if(sectionName == "pairs"){
                section = Pairs;
            }
            else if(sectionName == "angles"){
                section = Angles;
            }
            else if(sectionName == "dihedrals"){
                section = Dihedrals;
            }
            else{
                section = Unknown;
            }
        }
        else if(section == Atoms){
            boost::algorithm::trim_left(line);
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            size_t index = boost::lexical_cast<size_t>(tokens[0]) - 1;

            topology->resize(index + 1);

            topology->setType(index, tokens[1]);

            chemkit::Real charge = boost::lexical_cast<chemkit::Real>(tokens[6]);
            topology->setCharge(index, charge);

            chemkit::Real mass = boost::lexical_cast<chemkit::Real>(tokens[7]);
            topology->setMass(index, mass);
        }
        else if(section == Bonds){
            boost::algorithm::trim_left(line);
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            size_t indexA = boost::lexical_cast<size_t>(tokens[0]);
            size_t indexB = boost::lexical_cast<size_t>(tokens[1]);

            topology->addBondedInteraction(indexA, indexB);
        }
        else if(section == Pairs){
            boost::algorithm::trim_left(line);
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            size_t indexA = boost::lexical_cast<size_t>(tokens[0]);
            size_t indexB = boost::lexical_cast<size_t>(tokens[1]);

            topology->addNonbondedInteraction(indexA, indexB);
        }
        else if(section == Angles){
            boost::algorithm::trim_left(line);
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            size_t indexA = boost::lexical_cast<size_t>(tokens[0]);
            size_t indexB = boost::lexical_cast<size_t>(tokens[1]);
            size_t indexC = boost::lexical_cast<size_t>(tokens[2]);

            topology->addAngleInteraction(indexA, indexB, indexC);
        }
        else if(section == Dihedrals){
            boost::algorithm::trim_left(line);
            std::vector<std::string> tokens;
            boost::split(tokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            size_t indexA = boost::lexical_cast<size_t>(tokens[0]);
            size_t indexB = boost::lexical_cast<size_t>(tokens[1]);
            size_t indexC = boost::lexical_cast<size_t>(tokens[2]);
            size_t indexD = boost::lexical_cast<size_t>(tokens[3]);

            topology->addTorsionInteraction(indexA, indexB, indexC, indexD);
        }
    }

    file->setTopology(topology);

    return true;
}
