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

#include "cubefileformat.h"

#include <QtCore>

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
        position *= BohrToAnstroms;

        // set position
        atom->setPosition(position);
    }

    file->addMolecule(molecule);

    return true;
}
