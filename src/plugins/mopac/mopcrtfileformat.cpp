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

#include "mopcrtfileformat.h"

#include <chemkit/molecule.h>

MopcrtFileFormat::MopcrtFileFormat()
    : chemkit::ChemicalFileFormat("mopcrt")
{
}

bool MopcrtFileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
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
    for(;;){
        line = iodev->readLine();
        QStringList lineItems = line.split(' ', QString::SkipEmptyParts);
        if(lineItems.size() < 7){
            break;
        }

        chemkit::Atom *atom = molecule->addAtom(lineItems[0].toStdString());
        if(!atom){
            continue;
        }

        atom->setPosition(lineItems[1].toDouble(),
                          lineItems[3].toDouble(),
                          lineItems[5].toDouble());
    }

    if(molecule->isEmpty()){
        delete molecule;
        return false;
    }

    file->addMolecule(molecule);

    return true;
}
