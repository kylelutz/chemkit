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

#include "xyzfileformat.h"

#include <iomanip>

#include <chemkit/atom.h>
#include <chemkit/element.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

XyzFileFormat::XyzFileFormat()
    : chemkit::MoleculeFileFormat("xyz")
{
}

bool XyzFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // atom count line
    int atomCount = 0;
    input >> atomCount;
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // comment line (unused)
    std::string commentLine;
    std::getline(input, commentLine);
    CHEMKIT_UNUSED(commentLine);

    // create molecule
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);

    // read atoms and coordinates
    for(int i = 0; i < atomCount; i++){
        std::string symbol;
        double x = 0;
        double y = 0;
        double z = 0;

        input >> symbol >> x >> y >> z;
        if(input.fail()){
            input.clear();
        }

        // add atom from symbol or atomic number
        chemkit::Atom *atom = 0;
        if(symbol.empty()){
            continue;
        }
        else if(isdigit(symbol.at(0))){
            int atomicNumber = boost::lexical_cast<int>(symbol);
            atom = molecule->addAtom(atomicNumber);
        }
        else{
            atom = molecule->addAtom(symbol);
        }

        // set atom position
        if(atom){
            atom->setPosition(x, y, z);
        }
    }

    // add molecule to file
    file->addMolecule(molecule);

    return true;
}

bool XyzFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    boost::shared_ptr<chemkit::Molecule> molecule = file->molecule();
    if(!molecule){
        setErrorString("No molecule in file.");
        return false;
    }

    // atom count line
    output << molecule->atomCount() << "\n";

    // comment line
    output << "\n";

    // atoms and coordinates
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        output << std::showpoint
               << std::setw(3) << atom->symbol()
               << std::setw(15) << std::setprecision(5) << atom->x()
               << std::setw(15) << std::setprecision(5) << atom->y()
               << std::setw(15) << std::setprecision(5) << atom->z()
               << "\n";
    }

    return true;
}
