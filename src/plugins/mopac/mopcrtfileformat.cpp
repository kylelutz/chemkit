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

#include "mopcrtfileformat.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

MopcrtFileFormat::MopcrtFileFormat()
    : chemkit::MoleculeFileFormat("mopcrt")
{
}

bool MopcrtFileFormat::read(QIODevice *iodev, chemkit::MoleculeFile *file)
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
