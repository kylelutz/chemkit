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

#include "txyzfileformat.h"

#include <QtCore>

#include <chemkit/atom.h>
#include <chemkit/element.h>
#include <chemkit/moleculefile.h>

TxyzFileFormat::TxyzFileFormat()
    : chemkit::MoleculeFileFormat("txyz")
{
}

TxyzFileFormat::~TxyzFileFormat()
{
}

bool TxyzFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
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

    boost::shared_ptr<chemkit::Molecule> molecule(new chemkit::Molecule);
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
            QByteArray symbol = line[1].toAscii();
            atomicNumber = chemkit::Element(symbol.constData()).atomicNumber();
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

bool TxyzFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
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

bool TxyzFileFormat::write(const chemkit::MoleculeFile *file, QIODevice *iodev)
{
    iodev->setTextModeEnabled(true);

    if(file->isEmpty()){
        setErrorString("File is empty.");
        return false;
    }

    const boost::shared_ptr<chemkit::Molecule> &molecule = file->molecule();

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

        // list of atom indices for each bonded neighbor
        QList<int> neighborIndices;
        foreach(const chemkit::Atom *neighbor, atom->neighbors()){
            neighborIndices.append(neighbor->index());
        }
        qSort(neighborIndices);

        // write neighbor indices
        foreach(int neighborIndex, neighborIndices){
            iodev->write(QString("%1").arg(neighborIndex+1, 6).toAscii());
        }

        iodev->write("\n");

        index++;
    }

    return true;
}
