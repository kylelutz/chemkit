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

#include "fhzfileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/internalcoordinates.h>

FhzFileFormat::FhzFileFormat()
    : chemkit::ChemicalFileFormat("fhz")
{
}

FhzFileFormat::~FhzFileFormat()
{
}

bool FhzFileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
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

    chemkit::Molecule *molecule = new chemkit::Molecule;
    chemkit::InternalCoordinates coordinates(atomCount);

    for(int i = 0; i < atomCount; i++){
        // read atom line
        line = iodev->readLine();
        lineItems = line.split(' ', QString::SkipEmptyParts);
        if(lineItems.size() < 1){
            break;
        }

        // create atom
        chemkit::Atom *atom = molecule->addAtom(lineItems[0]);
        if(!atom){
            continue;
        }

        // read coordinates and connections
        int a = 0;
        int b = 0;
        int c = 0;
        chemkit::Float r = 0;
        chemkit::Float theta = 0;
        chemkit::Float phi = 0;

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
        coordinates.setCoordinates(i, r, theta, phi);
        coordinates.setConnections(i, a - 1, b - 1, c - 1);
    }

    molecule->setCoordinates(&coordinates);

    file->addMolecule(molecule);

    return true;
}
