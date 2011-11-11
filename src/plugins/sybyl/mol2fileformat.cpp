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

#include "mol2fileformat.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/moleculefile.h>

Mol2FileFormat::Mol2FileFormat()
    : chemkit::io::MoleculeFileFormat("mol2")
{
}

Mol2FileFormat::~Mol2FileFormat()
{
}

bool Mol2FileFormat::read(std::istream &input, chemkit::io::MoleculeFile *file)
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

bool Mol2FileFormat::read(QIODevice *iodev, chemkit::io::MoleculeFile *file)
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
                molecule->setName(name.toStdString());
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

                    chemkit::Atom *atom = molecule->addAtom(symbol.toStdString());
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

bool Mol2FileFormat::write(const chemkit::io::MoleculeFile *file, std::ostream &output)
{
    QBuffer buffer;
    buffer.open(QBuffer::WriteOnly);
    bool ok = write(file, &buffer);
    if(!ok){
        return false;
    }

    output.write(buffer.data().constData(), buffer.size());
    return true;
}

bool Mol2FileFormat::write(const chemkit::io::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    foreach(chemkit::Molecule *molecule, file->molecules()){
        iodev->write("@<TRIPOS>MOLECULE\n");
        iodev->write((molecule->name() + "\n").c_str());
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
                                                              atom->symbol().c_str(),
                                                              atomNumber,
                                                              atom->x(),
                                                              atom->y(),
                                                              atom->z(),
                                                              atom->symbol().c_str(),
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
