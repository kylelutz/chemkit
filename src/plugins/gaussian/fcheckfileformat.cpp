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

#include "fcheckfileformat.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/constants.h>

FCheckFileFormat::FCheckFileFormat()
    : chemkit::MoleculeFileFormat("fchk")
{
}

FCheckFileFormat::~FCheckFileFormat()
{
}

bool FCheckFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // find the number of atoms:
    int atomCount;
    std::string currentLine;
    while (std::getline(input, currentLine))
    {
        if (currentLine.find("Number of atoms") != std::string::npos)
        {
            std::istringstream iss(currentLine.substr(49));
            iss >> atomCount;
            break;
        }
    }

    // find the atomic numbers:
    std::vector<int> atomicNumbers(atomCount);
    while (std::getline(input, currentLine))
    {
        if (currentLine.find("Atomic numbers") != std::string::npos)
        {
            int remainder = atomCount % 6;

            for (int i = 0; i < atomCount - remainder; i += 6)
                input >> atomicNumbers[i + 0]
                      >> atomicNumbers[i + 1]
                      >> atomicNumbers[i + 2]
                      >> atomicNumbers[i + 3]
                      >> atomicNumbers[i + 4]
                      >> atomicNumbers[i + 5];

            for (int i = atomCount - remainder; i < atomCount; i++)
                input >> atomicNumbers[i];
            break;
        }
    }

    // find the positions:
    std::vector<double> positions(3 * atomCount);
    while (std::getline(input, currentLine))
    {
        if (currentLine.find("Current cartesian coordinates") != std::string::npos)
        {
            int remainder = (3 * atomCount) % 5;

            for (int i = 0; i < (3 * atomCount) - remainder; i += 5)
                input >> positions[i + 0]
                      >> positions[i + 1]
                      >> positions[i + 2]
                      >> positions[i + 3]
                      >> positions[i + 4];

            for (int i = (3 * atomCount) - remainder; i < (3 * atomCount); i++)
                input >> positions[i];
            break;
        }
    }

    // create the new molecule:
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    for (int i = 0; i < atomCount; i++)
    {
        double x = positions[3 * i + 0] * chemkit::constants::BohrToAngstroms;
        double y = positions[3 * i + 1] * chemkit::constants::BohrToAngstroms;
        double z = positions[3 * i + 2] * chemkit::constants::BohrToAngstroms;

        chemkit::Atom *atom = molecule->addAtom(atomicNumbers[i]);
        chemkit::Point3 position(x, y, z);
        atom->setPosition(position);
    }
    file->addMolecule(molecule);

    return true;
}
