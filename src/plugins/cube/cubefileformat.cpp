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

#include "cubefileformat.h"

#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/constants.h>
#include <chemkit/moleculefile.h>

CubeFileFormat::CubeFileFormat()
    : chemkit::MoleculeFileFormat("cube")
{
}

CubeFileFormat::~CubeFileFormat()
{
}

bool CubeFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    // title line
    std::string titleLine;
    std::getline(input, titleLine);
    std::vector<std::string> titleLineItems;
    boost::split(titleLineItems, titleLine, boost::is_any_of("\t "));
    if(titleLineItems.size() > 0){
        molecule->setName(titleLineItems[0]);
    }

    // comment line
    std::string commentLine;
    std::getline(input, commentLine);

    // atom count line
    int atomCount = 0;
    input >> atomCount;
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    atomCount = std::abs(atomCount);

    // voxel count and axii lines
    std::string line;
    std::getline(input, line);
    std::getline(input, line);
    std::getline(input, line);

    // atom lines
    for(int i = 0; i < atomCount; i++){
        // read atomic number
        int atomicNumber = 0;
        input >> atomicNumber;

        // add atom
        chemkit::Atom *atom = molecule->addAtom(atomicNumber);
        if(!atom){
            continue;
        }

        // read unused float value
        float value;
        input >> value;

        // read position
        chemkit::Real x = 0;
        chemkit::Real y = 0;
        chemkit::Real z = 0;
        input >> x >> y >> z;

        chemkit::Point3 position(x, y, z);

        // scale from bohr units to angstroms
        position *= chemkit::constants::BohrToAnstroms;

        // set position
        atom->setPosition(position);

        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    file->addMolecule(molecule);

    return true;
}
