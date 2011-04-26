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
