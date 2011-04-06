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

#include "txyzfileformat.h"

#include <chemkit/element.h>

TxyzFileFormat::TxyzFileFormat()
    : chemkit::MoleculeFileFormat("txyz")
{
}

TxyzFileFormat::~TxyzFileFormat()
{
}

bool TxyzFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
{
    iodev->setTextModeEnabled(true);

    bool ok;

    // first line contains atom count and molecule name
    QString firstLine = iodev->readLine();
    QStringList firstLineItems = firstLine.split(' ', QString::SkipEmptyParts);
    if(firstLineItems.size() < 1){
        setErrorString("First line of TXYZ file should contain number of atoms.");
        return false;
    }

    int atomCount = firstLineItems[0].toInt(&ok);
    if(!ok){
        setErrorString("First line of TXYZ file should contain number of atoms.");
        return false;
    }

    chemkit::Molecule *molecule = new chemkit::Molecule();
    QVector<QList<int> > bondLists(atomCount);

    for(int i = 0; i < atomCount; i++){
        QStringList line = QString(iodev->readLine()).simplified().split(' ', QString::SkipEmptyParts);
        if(line.size() < 5){
            // line too short
            continue;
        }

        int atomicNumber;

        // if first character of line is a number then it is the atomic number
        if(line[1][0].isNumber()){
            atomicNumber = line[0].toInt();
        }
        // else we interpret it as an atomic symbol
        else{
            atomicNumber = chemkit::Element(line[1].toStdString()).atomicNumber();
        }

        // add the atom if we have a valid atomic number
        if(chemkit::Element::isValidAtomicNumber(atomicNumber)){
            chemkit::Atom *atom = molecule->addAtom(atomicNumber);
            atom->setPosition(line[2].toDouble(), line[3].toDouble(), line[4].toDouble());
        }

        // read bonds
        for(int j = 6; j < line.size(); j++){
            // index of the other atom in the bond
            int otherAtom = line[j].toInt(&ok);
            if(ok)
                bondLists[i].append(otherAtom);
        }
    }

    // make bonds
    for(int i = 0; i < atomCount; i++){
        QList<int> bondedAtoms = bondLists[i];

        foreach(int neighbor, bondedAtoms){
            molecule->addBond(i, neighbor - 1);
        }
    }

    file->addMolecule(molecule);

    return true;
}

bool TxyzFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    if(file->isEmpty()){
        setErrorString("File is empty.");
        return false;
    }

    const chemkit::Molecule *molecule = file->molecule();

    // write atom count and molecule name
    iodev->write(QString("%1 %2\n").arg(molecule->atomCount())
                                   .arg(molecule->name().c_str())
                                   .toAscii());

    int index = 0;
    foreach(const chemkit::Atom *atom, molecule->atoms()){
        // write atom line: index, symbol, x, y, z, 0
        iodev->write(QString("%1%2%3%4%5%6").arg(index+1, 6)
                                            .arg(atom->symbol().c_str(), 3)
                                            .arg(atom->x(), 10)
                                            .arg(atom->y(), 10)
                                            .arg(atom->z(), 10)
                                            .arg("0", 6)
                                            .toAscii());

        // list of atom indicies for each bonded neighbor
        QList<int> neighborIndicies;
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            neighborIndicies.append(neighbor->index());
        }
        qSort(neighborIndicies);

        // write neighbor indicies
        foreach(int neighborIndex, neighborIndicies){
            iodev->write(QString("%1").arg(neighborIndex+1, 6).toAscii());
        }

        iodev->write("\n");

        index++;
    }

    return true;
}
