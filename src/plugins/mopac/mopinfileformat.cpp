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

#include "mopinfileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/coordinateset.h>
#include <chemkit/internalcoordinates.h>

MopinFileFormat::MopinFileFormat()
    : chemkit::MoleculeFileFormat("mopin")
{
}

bool MopinFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // create molecule
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    // keyword line
    std::string keyword;
    std::getline(input, keyword);

    // title line
    std::string title;
    std::getline(input, title);
    if(!title.empty()){
        molecule->setName(title);
    }

    // blank line
    std::string line;
    std::getline(input, line);

    // read atoms
    std::vector<chemkit::Real> coordinateValues;
    std::vector<size_t> connectionValues;
    for(;;){
        std::getline(input, line);
        boost::trim(line);
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
        if(tokens.size() < 9){
            break;
        }

        chemkit::Atom *atom = molecule->addAtom(tokens[0]);
        if(!atom){
            continue;
        }

        coordinateValues.push_back(boost::lexical_cast<chemkit::Real>(tokens[1]));
        coordinateValues.push_back(boost::lexical_cast<chemkit::Real>(tokens[3]));
        coordinateValues.push_back(boost::lexical_cast<chemkit::Real>(tokens[5]));

        connectionValues.push_back(boost::lexical_cast<size_t>(tokens[7]));
        connectionValues.push_back(boost::lexical_cast<size_t>(tokens[8]));
        connectionValues.push_back(boost::lexical_cast<size_t>(tokens[9]));
    }

    chemkit::InternalCoordinates *coordinates = new chemkit::InternalCoordinates(molecule->size());

    for(size_t i = 0; i < molecule->size(); i++){
        coordinates->setCoordinates(i, coordinateValues[i*3+0],
                                       coordinateValues[i*3+1],
                                       coordinateValues[i*3+2]);

        coordinates->setConnections(i, connectionValues[i*3+0] - 1,
                                       connectionValues[i*3+1] - 1,
                                       connectionValues[i*3+2] - 1);
    }

    // set molecule coordinates
    molecule->addCoordinateSet(coordinates);

    // add molecule to the file
    file->addMolecule(molecule);

    return true;
}
