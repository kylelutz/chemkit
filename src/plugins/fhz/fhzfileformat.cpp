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

#include "fhzfileformat.h"

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>
#include <chemkit/internalcoordinates.h>

FhzFileFormat::FhzFileFormat()
    : chemkit::MoleculeFileFormat("fhz")
{
}

FhzFileFormat::~FhzFileFormat()
{
}

bool FhzFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
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

bool FhzFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    // title line
    QString line = iodev->readLine();

    // atom count line
    line = iodev->readLine();
    QStringList lineItems = line.split(' ', QString::SkipEmptyParts);
    if(lineItems.size() < 1){
        setErrorString("Failed to read atom count.");
        return false;
    }
    int atomCount = lineItems[0].toInt();

    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
    chemkit::InternalCoordinates *coordinates = new chemkit::InternalCoordinates(atomCount);

    for(int i = 0; i < atomCount; i++){
        // read atom line
        line = iodev->readLine();
        lineItems = line.split(' ', QString::SkipEmptyParts);
        if(lineItems.size() < 1){
            break;
        }

        // create atom
        QByteArray symbol = lineItems[0].toAscii();
        chemkit::Atom *atom = molecule->addAtom(symbol.constData());
        if(!atom){
            continue;
        }

        // read coordinates and connections
        int a = 0;
        int b = 0;
        int c = 0;
        chemkit::Real r = 0;
        chemkit::Real theta = 0;
        chemkit::Real phi = 0;

        switch(i){
            default:
                c = lineItems[5].toInt();
                phi = lineItems[6].toDouble();
            case 2:
                b = lineItems[3].toInt();
                theta = lineItems[4].toDouble();
            case 1:
                a = lineItems[1].toInt();
                r = lineItems[2].toDouble();
            case 0:
                break;
        }

        // set coordinates
        coordinates->setCoordinates(i, r, theta, phi);
        coordinates->setConnections(i, a - 1, b - 1, c - 1);
    }

    molecule->addCoordinateSet(coordinates);

    file->addMolecule(molecule);

    return true;
}
