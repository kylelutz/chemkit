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

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

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
    QByteArray data;
    while(!input.eof()){
        data += input.get();
    }
    data.chop(1);

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);
    return read(&buffer, file);
}

bool MdlFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    if(name() == "mol" || name() == "mdl"){
        return readMolFile(iodev, file);
    }
    else if(name() == "sdf" || name() == "sd"){
        return readSdfFile(iodev, file);
    }
    else{
        return false;
    }
}

bool MdlFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    QBuffer buffer;
    buffer.open(QBuffer::WriteOnly);
    bool ok = write(file, &buffer);
    if(!ok){
        return false;
    }

    output.write(buffer.data().constData(), buffer.size());
    return true;
}

bool MdlFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    if(file->isEmpty()){
        return false;
    }

    iodev->setTextModeEnabled(true);

    if(name() == "mol" || name() == "mdl"){
        writeMolFile(file->molecule(), iodev);
    }
    else if(name() == "sdf" || name() == "sd"){
        writeSdfFile(file, iodev);
    }
    else{
        return false;
    }

    return true;
}

// --- Internal Methods ---------------------------------------------------- //
bool MdlFileFormat::readMolFile(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    QString title = iodev->readLine().trimmed(); // title line
    iodev->readLine(); // creator line
    iodev->readLine(); // comment line
    if(iodev->atEnd()){
        setErrorString("File is empty");
        return false;
    }

    // read counts line
    char countsLine[41];
    int countsLineLength = iodev->readLine(countsLine, sizeof(countsLine));
    if(countsLineLength == -1){
        setErrorString("Failed to read counts line.");
        return false;
    }

    int atomCount = 0;
    int bondCount = 0;
    sscanf(countsLine, "%3d%3d", &atomCount, &bondCount);

    // create molecule
    chemkit::Molecule *molecule = new chemkit::Molecule;
    if(!title.isEmpty()){
        molecule->setName(title.toStdString());
    }

    // read atoms
    readAtomBlock(iodev, molecule, atomCount);

    // read bonds
    readBondBlock(iodev, molecule, bondCount);

    // read properties
    readPropertyBlock(iodev, molecule);

    // add molecule to the file
    file->addMolecule(molecule);

    return true;
}

bool MdlFileFormat::readSdfFile(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    while(!iodev->atEnd()){
        // read molecule
        readMolFile(iodev, file);

        // read data block
        chemkit::Molecule *molecule = file->molecules().back();
        readDataBlock(iodev, molecule, file);
    }

    // return false if we failed to read any molecules
    if(file->moleculeCount() == 0){
        return false;
    }

    return true;
}

bool MdlFileFormat::readAtomBlock(QIODevice *iodev, chemkit::Molecule *molecule, int atomCount)
{
    for(int i = 0; i < atomCount; i++){
        char line[72];

        int lineLength = iodev->readLine(line, 72);
        if(lineLength == -1){
            return false;
        }
        if(lineLength < 33){
            // line too short
            continue;
        }

        double x, y, z;
        char symbol[3];
        sscanf(line, "%10lf%10lf%10lf%3s", &x, &y, &z, symbol);

        chemkit::Atom *atom = molecule->addAtom(symbol);
        if(!atom){
            // invalid atomic symbol
            return false;
        }

        atom->setPosition(x, y, z);
    }

    return true;
}

bool MdlFileFormat::readBondBlock(QIODevice *iodev, chemkit::Molecule *molecule, int bondCount)
{
    for(int i = 0; i < bondCount; i++){
        char line[24];
        int lineLength = iodev->readLine(line, 24);
        if(lineLength == -1){
            return false;
        }
        else if(lineLength < 9){
            // line too short
            return false;
        }

        int firstAtomIndex;
        int secondAtomIndex;
        int bondOrder;
        sscanf(line, "%3d%3d%3d", &firstAtomIndex, &secondAtomIndex, &bondOrder);

        chemkit::Atom *firstAtom = molecule->atom(firstAtomIndex - 1);
        chemkit::Atom *secondAtom = molecule->atom(secondAtomIndex - 1);
        if(firstAtom && secondAtom){
            chemkit::Bond *bond = molecule->addBond(firstAtom, secondAtom);

            // aromatic bond
            if(bondOrder == 4)
                bondOrder = 1;

            bond->setOrder(bondOrder);
        }
    }

    return true;
}

bool MdlFileFormat::readPropertyBlock(QIODevice *iodev, chemkit::Molecule *molecule)
{
    CHEMKIT_UNUSED(molecule);

    while(!iodev->atEnd()){
        QString line = iodev->readLine();

        if(line.startsWith("M  END")){
            return true;
        }
    }

    return false;
}

bool MdlFileFormat::readDataBlock(QIODevice *iodev, chemkit::Molecule *molecule, chemkit::MoleculeFile *file)
{
    QString dataName;
    QString dataValue;

    bool readingValue = false;

    while(!iodev->atEnd()){
        QString line = iodev->readLine().trimmed();

        if(line.startsWith("$$$$")){
            return true;
        }
        else if(line.startsWith("> <")){
            dataName = line.mid(3, line.length() - 4);
            readingValue = true;
        }
        else if(readingValue && line.isEmpty()){
            molecule->setData(dataName.toStdString(), dataValue.toStdString());
            dataValue.clear();
        }
        else if(readingValue){
            if(!dataValue.isEmpty())
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

void MdlFileFormat::writeMolFile(const chemkit::Molecule *molecule, QIODevice *iodev)
{
    // name, creator, and comment lines
    iodev->write(molecule->name().c_str());
    iodev->write("\n\n\n");

    // counts line
    char countsLine[41];
    sprintf(countsLine, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n", molecule->atomCount(),
                                                                     molecule->bondCount());
    iodev->write(countsLine, sizeof(countsLine) - 1);

    // atoms
    writeAtomBlock(molecule, iodev);

    // bonds
    writeBondBlock(molecule, iodev);

    // properties
    iodev->write("M  END\n");
}

void MdlFileFormat::writeSdfFile(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    foreach(const chemkit::Molecule *molecule, file->molecules()){
        writeMolFile(molecule, iodev);
        iodev->write("$$$$\n");
    }
}

void MdlFileFormat::writeAtomBlock(const chemkit::Molecule *molecule, QIODevice *iodev)
{
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        char line[50];
        sprintf(line, "%10.4f%10.4f%10.4f %3s 0  0  0  0  0\n", atom->x(),
                                                                atom->y(),
                                                                atom->z(),
                                                                atom->symbol().c_str());
        iodev->write(line, sizeof(line) - 1);
    }
}

void MdlFileFormat::writeBondBlock(const chemkit::Molecule *molecule, QIODevice *iodev)
{
    foreach(const chemkit::Bond *bond, molecule->bonds()){
        char line[23];
        sprintf(line, "%3d%3d%3d  0  0  0  0\n", bond->atom1()->index() + 1,
                                                 bond->atom2()->index() + 1,
                                                 bond->order());
        iodev->write(line, sizeof(line) - 1);
    }
}
