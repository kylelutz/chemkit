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

#include "mol2fileformat.h"

Mol2FileFormat::Mol2FileFormat()
    : chemkit::ChemicalFileFormat("mol2")
{
}

Mol2FileFormat::~Mol2FileFormat()
{
}

bool Mol2FileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
{
    iodev->setTextModeEnabled(true);

    chemkit::Molecule *molecule = 0;

    QHash<int, chemkit::Atom *> atom_ids;

    int atomCount = 0;
    int bondCount = 0;

    while(!iodev->atEnd()){
        QString line = iodev->readLine();

        if(line.startsWith("@<TRIPOS>MOLECULE")){
            if(molecule){
                file->addMolecule(molecule);
            }

            QString name = iodev->readLine().simplified();

            QStringList countsLine = QString(iodev->readLine()).split(' ', QString::SkipEmptyParts);
            if(countsLine.size() < 2){
                qDebug() << "chemkit: Mol2FileFormat: Counts line too short.";
                continue;
            }

            atomCount = countsLine[0].toInt();
            bondCount = countsLine[1].toInt();

            atom_ids.clear();
            molecule = new chemkit::Molecule();

            if(!name.isEmpty())
                molecule->setName(name);
        }
        else if(!molecule){
            continue;
        }
        else if(line.startsWith("@<TRIPOS>")){
            line = line.trimmed();
            QStringRef type = line.midRef(9);

            if(type == "ATOM"){
                while(atomCount--){
                    QStringList atomLine = QString(iodev->readLine()).split(' ', QString::SkipEmptyParts);
                    if(atomLine.size() < 6){
                        qDebug() << "chemkit: Mol2FileFormat: Atom line too short.";
                        return false;
                    }

                    QString symbol = QString(atomLine[5]).split('.', QString::SkipEmptyParts)[0];
                    if(symbol.length() > 1){
                        symbol = symbol[0] + symbol.mid(1).toLower();
                    }

                    chemkit::Atom *atom = molecule->addAtom(symbol);
                    if(!atom){
                        qDebug() << "chemkit: Mol2FileFormat: Invalid atom symbol: " << symbol;
                        continue;
                    }

                    atom->setPosition(atomLine[2].toDouble(),
                                      atomLine[3].toDouble(),
                                      atomLine[4].toDouble());

                    if(atomLine.size() >= 9){
                        atom->setPartialCharge(atomLine[8].toDouble());
                    }

                    atom_ids[atomLine[0].toInt()] = atom;
                }
            }
            else if(type == "BOND"){
                while(bondCount--){
                    QStringList bondLine = QString(iodev->readLine()).split(' ', QString::SkipEmptyParts);
                    if(bondLine.size() < 4){
                        qDebug() << "chemkit: Mol2FileFormat: Bond line too short.";
                        continue;
                    }

                    chemkit::Atom *a1 = atom_ids.value(bondLine[1].toInt(), 0);
                    chemkit::Atom *a2 = atom_ids.value(bondLine[2].toInt(), 0);

                    if(!a1 || !a2){
                        continue;
                    }

                    int order;
                    QString bondOrder = bondLine[3];

                    bool is_int;
                    order = bondOrder.toInt(&is_int);
                    if(!is_int){
                        // aromatic bond
                        if(bondOrder == "ar")
                            order = 1;
                        // amide bond
                        else if(bondOrder == "am")
                            order = 1;
                        // not connected
                        else if(bondOrder == "nc")
                            order = 0;
                        // default = single bond
                        else
                            order = 1;
                    }

                    if(order > 0){
                        molecule->addBond(a1, a2, order);
                    }
                }
            }
        }
    }

    if(molecule){
        file->addMolecule(molecule);
    }

    return true;
}

bool Mol2FileFormat::write(const chemkit::ChemicalFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    foreach(chemkit::Molecule *molecule, file->molecules()){
        iodev->write("@<TRIPOS>MOLECULE\n");
        iodev->write(molecule->name().toAscii() + "\n");
        QString countsLine;
        countsLine.sprintf("%4u%4u%3u%3u%3u\n", molecule->atomCount(), molecule->bondCount(), 0, 0, 0);
        iodev->write(countsLine.toAscii());
        iodev->write("SMALL\n");
        iodev->write("GASTEIGER\n");
        iodev->write("\n");
        iodev->write("\n");

        iodev->write("@<TRIPOS>ATOM\n");
        int atomNumber = 1;
        foreach(chemkit::Atom *atom, molecule->atoms()){
            QString line;
            line.sprintf("%9u %s%u %g %g %g %s %u <%u> %g\n", atomNumber,
                                                              qPrintable(atom->symbol()),
                                                              atomNumber,
                                                              atom->x(),
                                                              atom->y(),
                                                              atom->z(),
                                                              qPrintable(atom->symbol()),
                                                              1,
                                                              1,
                                                              atom->partialCharge());
            iodev->write(line.toAscii());
            atomNumber++;
        }

        iodev->write("@<TRIPOS>BOND\n");
        int bondNumber = 1;
        foreach(chemkit::Bond *bond, molecule->bonds()){
            QString line;
            line.sprintf("%6u%6u%6u%6u\n", bondNumber,
                                           bond->atom1()->index()+1,
                                           bond->atom2()->index()+1,
                                           bond->order());
            iodev->write(line.toAscii());
            bondNumber++;
        }
    }

    return true;
}
