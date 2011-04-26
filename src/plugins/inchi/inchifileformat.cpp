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

#include "inchifileformat.h"

#include <chemkit/lineformat.h>

// --- Construction and Destruction ---------------------------------------- //
InchiFileFormat::InchiFileFormat()
    : chemkit::MoleculeFileFormat("inchi")
{
}

InchiFileFormat::~InchiFileFormat()
{
}

// --- Input/Output -------------------------------------------------------- //
bool InchiFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
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

        std::string inchi = splitLine[0].toStdString();
        chemkit::Molecule *molecule = inchiFormat->read(inchi);
        if(!molecule)
            continue;

        if(splitLine.size() >= 2){
            QString name = splitLine[1];
            molecule->setName(name.toStdString());
        }

        file->addMolecule(molecule);
    }

    delete inchiFormat;

    return true;
}

bool InchiFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    chemkit::LineFormat *inchiFormat = chemkit::LineFormat::create("inchi");
    if(!inchiFormat){
        setErrorString("InChI line format not supported.");
        return false;
    }

    foreach(const chemkit::Molecule *molecule, file->molecules()){
        std::string inchi = inchiFormat->write(molecule);
        iodev->write(inchi.c_str());

        if(!molecule->name().empty()){
            iodev->write(" ");
            iodev->write(molecule->name().c_str());
        }

        iodev->write("\n");
    }

    delete inchiFormat;

    return true;
}
