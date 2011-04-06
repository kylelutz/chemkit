/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include "smifileformat.h"

// --- Construction and Destruction ---------------------------------------- //
SmiFileFormat::SmiFileFormat()
    : chemkit::MoleculeFileFormat("smi")
{

}

SmiFileFormat::~SmiFileFormat()
{
}

// --- Input/Output -------------------------------------------------------- //
bool SmiFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    chemkit::LineFormat *smilesFormat = chemkit::LineFormat::create("smiles");
    if(!smilesFormat){
        setErrorString("SMILES line format not supported.");
        return false;
    }

    while(!iodev->atEnd()){
        QString line = iodev->readLine().simplified();

        QStringList splitLine = line.split(' ', QString::SkipEmptyParts);
        if(splitLine.isEmpty()){
            continue;
        }

        std::string smiles = splitLine[0].toStdString();
        chemkit::Molecule *molecule = smilesFormat->read(smiles);
        if(!molecule){
            qDebug() << "Error with smiles: " << smiles.c_str();
            qDebug() << smilesFormat->errorString().c_str();
            continue;
        }

        if(splitLine.size() >= 2){
            QString name = line.mid(smiles.length()).trimmed();
            molecule->setName(name.toStdString());
        }

        file->addMolecule(molecule);
    }

    delete smilesFormat;

    return true;
}

bool SmiFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    chemkit::LineFormat *smilesFormat = chemkit::LineFormat::create("smiles");
    if(!smilesFormat){
        setErrorString("SMILES line format not supported.");
        return false;
    }

    foreach(const chemkit::Molecule *molecule, file->molecules()){
        std::string smiles = smilesFormat->write(molecule);
        iodev->write(smiles.c_str());

        if(!molecule->name().empty()){
            iodev->write(" ");
            iodev->write(molecule->name().c_str());
        }

        iodev->write("\n");
    }

    delete smilesFormat;

    return true;
}
