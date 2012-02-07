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

// The MdlFileFormat class implements reading and writing of the mol,
// mdl, sd, and sdf chemical files.
//
// Specification: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp

#include "mdlfileformat.h"

#include <boost/algorithm/string.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

namespace {

int readNumber(const char *s, int length)
{
    // skip any leading white space
    while(isspace(*s) && --length){
        s++;
    }

    int number = 0;

    while(length-- && isdigit(*s)){
        number = number * 10 + (*s++ - '0');
    }

    return number;
}

} // end anonymous namespace

// --- Construction and Destruction ---------------------------------------- //
MdlFileFormat::MdlFileFormat(const std::string &name)
    : chemkit::MoleculeFileFormat(name)
{
}

MdlFileFormat::~MdlFileFormat()
{
}

// --- Input and Output ---------------------------------------------------- //
bool MdlFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    if(name() == "mol" || name() == "mdl"){
        return readMolFile(input, file);
    }
    else if(name() == "sdf" || name() == "sd"){
        return readSdfFile(input, file);
    }
    else{
        return false;
    }
}

bool MdlFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    if(file->isEmpty()){
        return false;
    }

    if(name() == "mol" || name() == "mdl"){
        writeMolFile(file->molecule().get(), output);
    }
    else if(name() == "sdf" || name() == "sd"){
        writeSdfFile(file, output);
    }
    else{
        return false;
    }

    return true;
}

// --- Internal Methods ---------------------------------------------------- //
bool MdlFileFormat::readMolFile(std::istream &input, chemkit::MoleculeFile *file)
{
    // title line
    std::string title;
    std::getline(input, title);

     // creator line
    std::string creator;
    std::getline(input, creator);

    // comment line
    std::string comment;
    std::getline(input, comment);

    if(input.eof()){
        setErrorString("File is empty");
        return false;
    }

    // read counts line
    std::string countsLine;
    std::getline(input, countsLine);

    int atomCount = readNumber(&countsLine[0], 3);
    int bondCount = readNumber(&countsLine[3], 3);

    // create molecule
    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    if(!title.empty()){
        molecule->setName(title);
    }

    // read atoms
    readAtomBlock(input, molecule.get(), atomCount);

    // read bonds
    readBondBlock(input, molecule.get(), bondCount);

    // read properties
    readPropertyBlock(input, molecule.get());

    file->addMolecule(molecule);

    return true;
}

bool MdlFileFormat::readSdfFile(std::istream &input, chemkit::MoleculeFile *file)
{
    while(!input.eof()){
        // read molecule
        readMolFile(input, file);

        // read data block
        chemkit::Molecule *molecule = file->molecules().back().get();
        readDataBlock(input, molecule);
    }

    // return false if we failed to read any molecules
    if(file->moleculeCount() == 0){
        return false;
    }

    return true;
}

bool MdlFileFormat::readAtomBlock(std::istream &input, chemkit::Molecule *molecule, int atomCount)
{
    for(int i = 0; i < atomCount; i++){
        std::string line;
        std::getline(input, line);

        if(line.size() < 33){
            // line too short
            continue;
        }

        double x, y, z;
        char symbol[3];
        sscanf(&line[0], "%10lf%10lf%10lf%3s", &x, &y, &z, symbol);

        chemkit::Atom *atom = molecule->addAtom(symbol);
        if(!atom->element().isValid()){
            if(strcmp(symbol, "D") == 0){
                atom->setIsotope(chemkit::Isotope(chemkit::Atom::Hydrogen, 2));
            }
            else if(strcmp(symbol, "T") == 0){
                atom->setIsotope(chemkit::Isotope(chemkit::Atom::Hydrogen, 3));
            }
        }

        atom->setPosition(x, y, z);
    }

    return true;
}

bool MdlFileFormat::readBondBlock(std::istream &input, chemkit::Molecule *molecule, int bondCount)
{
    for(int i = 0; i < bondCount; i++){
        std::string line;
        std::getline(input, line);
        if(line.size() < 9){
            // line too short
            return false;
        }

        int firstAtomIndex = readNumber(&line[0], 3);
        int secondAtomIndex = readNumber(&line[3], 3);

        chemkit::Atom *firstAtom = molecule->atom(firstAtomIndex - 1);
        chemkit::Atom *secondAtom = molecule->atom(secondAtomIndex - 1);
        if(firstAtom && secondAtom){
            chemkit::Bond *bond = molecule->addBond(firstAtom, secondAtom);

            int bondOrder = line[8] - '0';
            bond->setOrder(bondOrder);
        }
    }

    return true;
}

bool MdlFileFormat::readPropertyBlock(std::istream &input, chemkit::Molecule *molecule)
{
    CHEMKIT_UNUSED(molecule);

    while(!input.eof()){
        std::string line;
        std::getline(input, line);

        if(boost::starts_with(line, "M  END")){
            return true;
        }
    }

    return false;
}

bool MdlFileFormat::readDataBlock(std::istream &input, chemkit::Molecule *molecule)
{
    std::string dataName;
    std::string dataValue;

    bool readingValue = false;

    while(!input.eof()){
        std::string line;
        std::getline(input, line);
        boost::algorithm::trim(line);

        if(boost::starts_with(line, "$$$$")){
            return true;
        }
        else if(boost::starts_with(line, "> <")){
            dataName = line.substr(3, line.length() - 4);
            readingValue = true;
        }
        else if(readingValue && line.empty()){
            molecule->setData(dataName, dataValue);
            dataValue.clear();
        }
        else if(readingValue){
            if(!dataValue.empty())
                dataValue += "\n";
            dataValue += line;
        }
        else{
            dataValue = line;
            readingValue = false;
        }
    }

    return false;
}

void MdlFileFormat::writeMolFile(const chemkit::Molecule *molecule, std::ostream &output)
{
    // name, creator, and comment lines
    output << molecule->name() << "\n";
    output << "\n";
    output << "\n";

    // counts line
    char countsLine[41];
    sprintf(countsLine, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n", static_cast<int>(molecule->atomCount()),
                                                                     static_cast<int>(molecule->bondCount()));
    output.write(countsLine, sizeof(countsLine) - 1);

    // atoms
    writeAtomBlock(molecule, output);

    // bonds
    writeBondBlock(molecule, output);

    // properties
    output << "M  END\n";
}

void MdlFileFormat::writeSdfFile(const chemkit::MoleculeFile *file, std::ostream &output)
{
    foreach(const boost::shared_ptr<chemkit::Molecule> molecule, file->molecules()){
        writeMolFile(molecule.get(), output);
        output << "$$$$\n";
    }
}

void MdlFileFormat::writeAtomBlock(const chemkit::Molecule *molecule, std::ostream &output)
{
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        char line[50];
        sprintf(line, "%10.4f%10.4f%10.4f %3s 0  0  0  0  0\n", atom->x(),
                                                                atom->y(),
                                                                atom->z(),
                                                                atom->symbol().c_str());
        output.write(line, sizeof(line) - 1);
    }
}

void MdlFileFormat::writeBondBlock(const chemkit::Molecule *molecule, std::ostream &output)
{
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        char line[23];
        sprintf(line, "%3d%3d%3d  0  0  0  0\n", static_cast<int>(bond->atom1()->index()) + 1,
                                                 static_cast<int>(bond->atom2()->index()) + 1,
                                                 static_cast<int>(bond->order()));
        output.write(line, sizeof(line) - 1);
    }
}
