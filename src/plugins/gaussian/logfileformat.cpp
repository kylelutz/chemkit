/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** This file is Copyright (C) 2016 by Jan von Cosel
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

#include "logfileformat.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

#include <iostream>

LogFileFormat::LogFileFormat()
    : chemkit::MoleculeFileFormat("log")
{
}

LogFileFormat::~LogFileFormat()
{
}

bool LogFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // find the number of atoms:
    int atomCount;
    std::string currentLine;
    while (std::getline(input, currentLine))
    {
        if (currentLine.find("NAtoms=") != std::string::npos)
        {
            std::istringstream iss(currentLine.substr(8));
            iss >> atomCount;
            break;
        }
    }

    // read in every Standard Orientation in the file and add them as molecules:
    while (std::getline(input, currentLine))
    {
        if (currentLine.find("Standard orientation") != std::string::npos)
        {
            boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

            // skip the next 4 lines:
            for (int i = 0; i < 4; i++)
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            int dummy1, atomicNumber, dummy2;
            double x, y, z;

            for (int i = 0; i < atomCount; i++)
            {
                input >> dummy1 >> atomicNumber >> dummy2 >> x >> y >> z;

                chemkit::Atom *atom = molecule->addAtom(atomicNumber);
                chemkit::Point3 position(x, y, z);
                atom->setPosition(position);
            }
            file->addMolecule(molecule);
        }
    }

    return true;
}
