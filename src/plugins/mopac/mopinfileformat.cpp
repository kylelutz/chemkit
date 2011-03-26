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

#include "mopinfileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/internalcoordinates.h>

MopinFileFormat::MopinFileFormat()
    : chemkit::ChemicalFileFormat("mopin")
{
}

bool MopinFileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
{
    iodev->setTextModeEnabled(true);

    // create molecule
    chemkit::Molecule *molecule = new chemkit::Molecule;

    // keyword line
    QString line = iodev->readLine();

    // title line
    line = iodev->readLine();
    QString name = line.trimmed();
    molecule->setName(name.toStdString());

    // blank line
    iodev->readLine();

    // read atoms
    QList<chemkit::Float> coordinateValues;
    QList<int> connectionValues;
    for(;;){
        line = iodev->readLine();
        QStringList lineItems = line.split(' ', QString::SkipEmptyParts);
        if(lineItems.size() < 9){
            break;
        }

        chemkit::Atom *atom = molecule->addAtom(lineItems[0].toStdString());
        if(!atom){
            continue;
        }

        coordinateValues.append(lineItems[1].toDouble());
        coordinateValues.append(lineItems[3].toDouble());
        coordinateValues.append(lineItems[5].toDouble());

        connectionValues.append(lineItems[7].toInt());
        connectionValues.append(lineItems[8].toInt());
        connectionValues.append(lineItems[9].toInt());
    }

    chemkit::InternalCoordinates coordinates(molecule->size());

    for(int i = 0; i < molecule->size(); i++){
        coordinates.setCoordinates(i, coordinateValues[i*3+0],
                                      coordinateValues[i*3+1],
                                      coordinateValues[i*3+2]);

        coordinates.setConnections(i, connectionValues[i*3+0] - 1,
                                      connectionValues[i*3+1] - 1,
                                      connectionValues[i*3+2] - 1);
    }

    // set molecule coordinates
    molecule->setCoordinates(&coordinates);

    file->addMolecule(molecule);

    return true;
}
