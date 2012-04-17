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

#include "fhzfileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/internalcoordinates.h>

FhzFileFormat::FhzFileFormat()
    : chemkit::MoleculeFileFormat("fhz")
{
}

bool FhzFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // title line
    std::string title;
    std::getline(input, title);

    // atom count line
    std::string line;
    std::getline(input, line);
    boost::trim(line);
    std::vector<std::string> tokens;
    boost::split(tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
    if(tokens.size() < 1){
        setErrorString("Failed to read atom count.");
        return false;
    }
    size_t atomCount = boost::lexical_cast<size_t>(tokens[0]);

    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    chemkit::InternalCoordinates *coordinates = new chemkit::InternalCoordinates(atomCount);

    for(size_t i = 0; i < atomCount; i++){
        // read atom line
        std::getline(input, line);
        boost::trim(line);
        boost::split(tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
        if(tokens.size() < 1){
            break;
        }

        // create atom
        chemkit::Atom *atom = molecule->addAtom(tokens[0]);
        if(!atom){
            continue;
        }

        // read coordinates and connections
        size_t a = 0;
        size_t b = 0;
        size_t c = 0;
        chemkit::Real r = 0;
        chemkit::Real theta = 0;
        chemkit::Real phi = 0;

        switch(i){
            default:
                c = boost::lexical_cast<size_t>(tokens[5]);
                phi = boost::lexical_cast<chemkit::Real>(tokens[6]);
            case 2:
                b = boost::lexical_cast<size_t>(tokens[3]);
                theta = boost::lexical_cast<chemkit::Real>(tokens[4]);
            case 1:
                a = boost::lexical_cast<size_t>(tokens[1]);
                r = boost::lexical_cast<chemkit::Real>(tokens[2]);
            case 0:
                break;
        }

        // set coordinates
        coordinates->setCoordinates(i, r, theta, phi);
        coordinates->setConnections(i, a - 1, b - 1, c - 1);
    }

    molecule->addCoordinateSet(coordinates);

    file->addMolecule(molecule);

    return true;
}
