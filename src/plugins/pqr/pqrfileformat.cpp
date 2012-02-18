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

#include "pqrfileformat.h"

#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/element.h>
#include <chemkit/foreach.h>
#include <chemkit/moleculefile.h>

PqrFileFormat::PqrFileFormat()
    : chemkit::MoleculeFileFormat("pqr")
{
}

PqrFileFormat::~PqrFileFormat()
{
}

bool PqrFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    for(;;){
        std::string line;
        std::getline(input, line);

        if(line.empty()){
            break;
        }

        if(boost::starts_with(line, "ATOM")){
            // split string by whitespace
            std::vector<std::string> lineTokens;
            boost::split(lineTokens,
                         line,
                         boost::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            if(lineTokens.size() < 10){
                break;
            }

            // read symbol as first character of atom name
            char symbol = lineTokens[2][0];

            // add atom
            chemkit::Atom *atom = molecule->addAtom(chemkit::Element::fromSymbol(symbol));

            // set coordinates
            atom->setPosition(boost::lexical_cast<chemkit::Real>(lineTokens[5]),
                              boost::lexical_cast<chemkit::Real>(lineTokens[6]),
                              boost::lexical_cast<chemkit::Real>(lineTokens[7]));

            // set charge
            atom->setPartialCharge(boost::lexical_cast<chemkit::Real>(lineTokens[8]));
        }
    }

    file->addMolecule(molecule);

    return true;
}
