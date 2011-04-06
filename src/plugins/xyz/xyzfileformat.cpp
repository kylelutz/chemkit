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

#include "xyzfileformat.h"

#include <chemkit/element.h>

XyzFileFormat::XyzFileFormat()
    : chemkit::MoleculeFileFormat("xyz")
{
}

XyzFileFormat::~XyzFileFormat()
{
}

bool XyzFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    bool ok;
    int atoms_count = iodev->readLine().trimmed().toInt(&ok);
    if(!ok){
        setErrorString("First line of XYZ file should contain number of atoms.");
        return false;
    }

    QString comment_line = iodev->readLine();
    Q_UNUSED(comment_line);

    chemkit::Molecule *molecule = new chemkit::Molecule();

    for(int i = 0; i < atoms_count; i++){
        QStringList line = QString(iodev->readLine()).simplified().split(' ', QString::SkipEmptyParts);
        if(line.size() < 4){
            // line too short
            continue;
        }

        int atomicNumber;

        // if first character of line is a number then it is the atomic number
        if(line[0][0].isNumber()){
            atomicNumber = line[0].toInt();
        }
        // else we interpret it as an atomic symbol
        else{
            atomicNumber = chemkit::Element(line[0].toStdString()).atomicNumber();
        }

        // add the atom if we have a valid atomic number
        if(chemkit::Element::isValidAtomicNumber(atomicNumber)){
            chemkit::Atom *atom = molecule->addAtom(atomicNumber);
            atom->setPosition(line[1].toDouble(), line[2].toDouble(), line[3].toDouble());
        }
    }

    file->addMolecule(molecule);

    return true;
}

bool XyzFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    chemkit::Molecule *molecule = file->molecule();
    if(!molecule){
        setErrorString("No molecule in file.");
        return false;
    }

    iodev->write(QString::number(molecule->atomCount()).toAscii() + "\n"); // atom count line
    iodev->write("\n");                                                    // comment line

    foreach(chemkit::Atom *atom, molecule->atoms()){
        QString line = QString("%1 %2 %3 %4\n").arg(atom->symbol().c_str())
                                               .arg(QString::number(atom->x()))
                                               .arg(QString::number(atom->y()))
                                               .arg(QString::number(atom->z()));
        iodev->write(line.toAscii());
    }

    return true;
}
