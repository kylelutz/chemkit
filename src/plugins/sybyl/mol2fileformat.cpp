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

#include "mol2fileformat.h"

#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

#include "sybylatomtyper.h"

Mol2FileFormat::Mol2FileFormat()
    : chemkit::MoleculeFileFormat("mol2")
{
}

Mol2FileFormat::~Mol2FileFormat()
{
}

bool Mol2FileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    boost::shared_ptr<chemkit::Molecule> molecule;

    int atomCount = 0;
    int bondCount = 0;

    while(!input.eof()){
        std::string line;
        std::getline(input, line);

        if(boost::starts_with(line, "@<TRIPOS>MOLECULE")){
            if(molecule){
                file->addMolecule(molecule);
            }

            std::string name;
            std::getline(input, name);
            boost::trim(name);

            std::string countsLineString;
            std::getline(input, countsLineString);
            boost::trim_left(countsLineString);
            std::vector<std::string> countsLine;
            boost::split(countsLine,
                         countsLineString,
                         boost::is_any_of(" \t"),
                         boost::token_compress_on);
            if(countsLine.size() < 2){
                continue;
            }

            atomCount = boost::lexical_cast<int>(countsLine[0]);
            bondCount = boost::lexical_cast<int>(countsLine[1]);

            molecule = boost::make_shared<chemkit::Molecule>();

            if(!name.empty()){
                molecule->setName(name);
            }
        }
        else if(!molecule){
            continue;
        }
        else if(boost::starts_with(line, "@<TRIPOS>")){
            boost::trim(line);
            std::string type = line.substr(9);

            if(type == "ATOM"){
                while(atomCount--){
                    std::string atomLineString;
                    std::getline(input, atomLineString);
                    boost::trim_left(atomLineString);
                    std::vector<std::string> atomLine;
                    boost::split(atomLine,
                                 atomLineString,
                                 boost::is_any_of(" \t"),
                                 boost::token_compress_on);
                    if(atomLine.size() < 6){
                        return false;
                    }

                    chemkit::Element element;
                    std::string symbol = atomLine[5];
                    boost::to_lower(symbol);
                    symbol[0] = toupper(symbol[0]);
                    size_t dotPosition = symbol.find_first_of('.');
                    if(dotPosition == std::string::npos){
                        element = chemkit::Element::fromSymbol(symbol);
                    }
                    else{
                        element = chemkit::Element::fromSymbol(symbol.c_str(), dotPosition);
                    }

                    chemkit::Atom *atom = molecule->addAtom(element);
                    if(!atom){
                        continue;
                    }

                    atom->setPosition(boost::lexical_cast<chemkit::Real>(atomLine[2]),
                                      boost::lexical_cast<chemkit::Real>(atomLine[3]),
                                      boost::lexical_cast<chemkit::Real>(atomLine[4]));

                    if(atomLine.size() >= 9){
                        atom->setPartialCharge(boost::lexical_cast<chemkit::Real>(atomLine[8]));
                    }
                }
            }
            else if(type == "BOND"){
                while(bondCount--){
                    std::string bondLineString;
                    std::getline(input, bondLineString);
                    boost::trim_left(bondLineString);
                    std::vector<std::string> bondLine;
                    boost::split(bondLine,
                                 bondLineString,
                                 boost::is_any_of(" \t"),
                                 boost::token_compress_on);
                    if(bondLine.size() < 4){
                        continue;
                    }

                    chemkit::Atom *a1 = molecule->atom(boost::lexical_cast<int>(bondLine[1]) - 1);
                    chemkit::Atom *a2 = molecule->atom(boost::lexical_cast<int>(bondLine[2]) - 1);

                    int bondOrder;
                    std::string bondOrderString = bondLine[3];

                    try {
                        bondOrder = boost::lexical_cast<int>(bondOrderString);
                    }
                    catch(boost::bad_lexical_cast&){
                        // aromatic bond
                        if(bondOrderString == "ar")
                            bondOrder = 1;
                        // amide bond
                        else if(bondOrderString == "am")
                            bondOrder = 1;
                        // not connected
                        else if(bondOrderString == "nc")
                            bondOrder = 0;
                        // default = single bond
                        else
                            bondOrder = 1;
                    }

                    if(bondOrder > 0){
                        molecule->addBond(a1, a2, bondOrder);
                    }
                }
            }
        }
    }

    if(molecule){
        file->addMolecule(molecule);
    }

    return true;
}

bool Mol2FileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    char line[80];

    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, file->molecules()){
        // perceive sybyl atom types
        SybylAtomTyper atomTyper;
        atomTyper.setMolecule(molecule.get());

        output << "@<TRIPOS>MOLECULE\n";
        output << molecule->name() << "\n";
        sprintf(line,
                "%4u%4u%3u%3u%3u\n",
                static_cast<int>(molecule->atomCount()),
                static_cast<int>(molecule->bondCount()),
                0,
                0,
                0);
        output << line;
        output << "SMALL\n";
        output << "GASTEIGER\n";
        output << "\n";
        output << "\n";

        output << "@<TRIPOS>ATOM\n";
        int atomNumber = 1;
        foreach(chemkit::Atom *atom, molecule->atoms()){
            // get atom type from the atom typer
            std::string type = atomTyper.typeString(atom);

            // use the atom's symbol if no type assigned
            if(type.empty()){
                type = atom->symbol();
            }

            sprintf(line,
                    "%7u %2s %10.4g %10.4g %10.4g %6s %u  LIG1 %10.4g\n",
                    atomNumber,
                    atom->symbol().c_str(),
                    atom->x(),
                    atom->y(),
                    atom->z(),
                    type.c_str(),
                    1,
                    atom->partialCharge());

            output << line;
            atomNumber++;
        }

        output << "@<TRIPOS>BOND\n";
        int bondNumber = 1;
        foreach(chemkit::Bond *bond, molecule->bonds()){
            sprintf(line,
                    "%6u%6u%6u%6u\n",
                    bondNumber,
                    static_cast<int>(bond->atom1()->index()) + 1,
                    static_cast<int>(bond->atom2()->index()) + 1,
                    bond->order());

            output << line;
            bondNumber++;
        }
    }

    return true;
}
