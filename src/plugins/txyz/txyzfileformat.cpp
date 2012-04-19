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

#include "txyzfileformat.h"

#include <iomanip>
#include <algorithm>

#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/element.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

TxyzFileFormat::TxyzFileFormat()
    : chemkit::MoleculeFileFormat("txyz")
{
}

bool TxyzFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // first line contains atom count and molecule name
    std::string firstLine;
    std::getline(input, firstLine);
    boost::trim(firstLine);
    std::vector<std::string> firstLineTokens;
    boost::split(firstLineTokens,
                 firstLine,
                 boost::is_any_of(" "),
                 boost::token_compress_on);
    if(firstLineTokens.size() < 1){
        setErrorString("First line of TXYZ file should contain number of atoms.");
        return false;
    }

    size_t atomCount = 0;
    try {
        atomCount = boost::lexical_cast<size_t>(firstLineTokens[0]);
    }
    catch(boost::bad_lexical_cast &e) {
        setErrorString("First line of TXYZ file should contain number of atoms.");
        return false;
    }

    // create molecule
    boost::shared_ptr<chemkit::Molecule> molecule = boost::make_shared<chemkit::Molecule>();

    // set molecule name
    if(firstLineTokens.size() >= 2){
        molecule->setName(firstLineTokens[1]);
    }

    // reserve space for atoms
    molecule->setAtomCapacity(atomCount);

    std::vector<std::vector<size_t> > bondLists(atomCount);

    for(size_t i = 0; i < atomCount; i++){
        std::string line;
        std::getline(input, line);
        boost::trim(line);
        std::vector<std::string> lineTokens;
        boost::split(lineTokens,
                     line,
                     boost::is_any_of("\t "),
                     boost::token_compress_on);
        if(lineTokens.size() < 5){
            // line too short
            continue;
        }

        chemkit::Atom::AtomicNumberType atomicNumber;

        // if first character of line is a number then it is the atomic number
        if(isdigit(lineTokens[1][0])){
            atomicNumber = boost::lexical_cast<chemkit::Atom::AtomicNumberType>(lineTokens[1]);
        }
        // else we interpret it as an atomic symbol
        else{
            atomicNumber = chemkit::Element(lineTokens[1]).atomicNumber();
        }

        // add the atom if we have a valid atomic number
        chemkit::Atom *atom = molecule->addAtom(atomicNumber);
        if(atom){
            atom->setPosition(boost::lexical_cast<chemkit::Real>(lineTokens[2]),
                              boost::lexical_cast<chemkit::Real>(lineTokens[3]),
                              boost::lexical_cast<chemkit::Real>(lineTokens[4]));
        }

        // read bonds
        for(size_t j = 6; j < lineTokens.size(); j++){
            try {
                bondLists[i].push_back(boost::lexical_cast<size_t>(lineTokens[j]));
            }
            catch(boost::bad_lexical_cast &){
                continue;
            }
        }
    }

    // add bonds
    for(size_t i = 0; i < atomCount; i++){
        const std::vector<size_t> &neighbors = bondLists[i];

        foreach(size_t neighbor, neighbors){
            molecule->addBond(i, neighbor - 1);
        }
    }

    file->addMolecule(molecule);

    return true;
}

bool TxyzFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    if(file->isEmpty()){
        setErrorString("File is empty.");
        return false;
    }

    const boost::shared_ptr<chemkit::Molecule> &molecule = file->molecule();

    // write atom count and molecule name
    output << std::setw(6) << molecule->atomCount();
    if(!molecule->name().empty()){
        output << "   " << molecule->name();
    }
    output << "\n";

    foreach(const chemkit::Atom *atom, molecule->atoms()){
        // write atom line: index, symbol, x, y, z, 0
        output << std::setw(6) << atom->index() + 1 <<
                  std::setw(4) << atom->symbol()  <<
                  std::setw(12) << atom->x() <<
                  std::setw(12) << atom->y() <<
                  std::setw(12) << atom->z() <<
                  std::setw(6) << "0";

        // list of atom indices for each bonded neighbor
        std::vector<size_t> neighborIndices;
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            neighborIndices.push_back(neighbor->index());
        }
        std::sort(neighborIndices.begin(), neighborIndices.end());

        // write neighbor indices
        foreach(size_t neighborIndex, neighborIndices){
            output << std::setw(6) << neighborIndex + 1;
        }

        output << "\n";
    }

    return true;
}
