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

#include "inchifileformat.h"

#include <chemkit/lineformat.h>

// --- Construction and Destruction ---------------------------------------- //
InchiFileFormat::InchiFileFormat()
    : chemkit::ChemicalFileFormat("inchi")
{
}

InchiFileFormat::~InchiFileFormat()
{
}

// --- Input/Output -------------------------------------------------------- //
bool InchiFileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
{
    iodev->setTextModeEnabled(true);

    chemkit::LineFormat *inchiFormat = chemkit::LineFormat::create("inchi");
    if(!inchiFormat){
        setErrorString("InChI line format not supported.");
        return false;
    }

    while(!iodev->atEnd()){
        QString line = iodev->readLine().trimmed();

        QStringList splitLine = line.split(" ");
        if(splitLine.isEmpty()){
            continue;
        }

        QString inchi = splitLine[0];
        chemkit::Molecule *molecule = inchiFormat->read(inchi);
        if(!molecule)
            continue;

        if(splitLine.size() >= 2){
            QString name = splitLine[1];
            molecule->setName(name);
        }

        file->addMolecule(molecule);
    }

    delete inchiFormat;

    return true;
}

bool InchiFileFormat::write(const chemkit::ChemicalFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    chemkit::LineFormat *inchiFormat = chemkit::LineFormat::create("inchi");
    if(!inchiFormat){
        setErrorString("InChI line format not supported.");
        return false;
    }

    foreach(const chemkit::Molecule *molecule, file->molecules()){
        QString inchi = inchiFormat->write(molecule);
        iodev->write(inchi.toAscii());

        if(molecule->hasName()){
            iodev->write(" ");
            iodev->write(molecule->name().toAscii());
        }

        iodev->write("\n");
    }

    delete inchiFormat;

    return true;
}
