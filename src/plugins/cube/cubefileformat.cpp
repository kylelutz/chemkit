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

#include "cubefileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

// Conversion factor between bohr units and Angstroms.
const double BohrToAnstroms = 0.52918;

CubeFileFormat::CubeFileFormat()
    : chemkit::MoleculeFileFormat("cube")
{
}

CubeFileFormat::~CubeFileFormat()
{
}

bool CubeFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    chemkit::Molecule *molecule = new chemkit::Molecule;

    // title line
    QString line = iodev->readLine();
    QStringList lineItems = line.split(" ", QString::SkipEmptyParts);
    if(lineItems.size() > 0){
        molecule->setName(lineItems[0].toStdString());
    }

    // comment line
    iodev->readLine();

    // atom count line
    line = iodev->readLine();
    lineItems = line.split(" ", QString::SkipEmptyParts);
    if(lineItems.size() < 1){
        setErrorString("Atom count line too short.");
        delete molecule;
        return false;
    }

    int atomCount = qAbs(lineItems[0].toInt());

    // voxel count and axii lines
    iodev->readLine();
    iodev->readLine();
    iodev->readLine();

    // atom lines
    for(int i = 0; i < atomCount; i++){
        line = iodev->readLine();
        lineItems = line.split(" ", QString::SkipEmptyParts);
        if(lineItems.size() < 5){
            continue;
        }

        int atomicNumber = lineItems[0].toInt();
        chemkit::Atom *atom = molecule->addAtom(atomicNumber);
        if(!atom){
            continue;
        }

        // read position
        chemkit::Point3 position(lineItems[2].toDouble(),
                                lineItems[3].toDouble(),
                                lineItems[4].toDouble());

        // scale from bohr units to angstroms
        position.scale(BohrToAnstroms);

        // set position
        atom->setPosition(position);
    }

    file->addMolecule(molecule);

    return true;
}
