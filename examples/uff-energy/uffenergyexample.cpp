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

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>
#include <chemkit/moleculefile.h>

int main(int argc, char *argv[])
{
    QTextStream out(stdout);
    QTextStream err(stderr);

    if(argc < 2){
        err << "Usage: " << argv[0] << " FILENAME" << "\n";
        return -1;
    }

    std::string fileName = argv[1];

    chemkit::MoleculeFile file(fileName);
    bool ok = file.read();
    if(!ok){
        err << "Failed to read file: " << fileName.c_str() << "\n";
        return -1;
    }

    boost::shared_ptr<chemkit::Molecule> molecule = file.molecule();
    if(!molecule){
        err << "File contains no molecules.\n";
        return -1;
    }

    chemkit::ForceField *uff = chemkit::ForceField::create("uff");
    if(!uff){
        err << "UFF force field plugin not found.\n";
        return -1;
    }

    uff->setMolecule(molecule.get());
    uff->setup();

    if(!uff->isSetup()){
        err << "Failed to parameterize force field.\n";
        delete uff;
        return -1;
    }

    chemkit::Real energy = uff->energy();

    out << "Formula: " << molecule->formula().c_str() << "\n";
    out << "Energy: " << energy << " kcal/mol\n";

    delete uff;

    return 0;
}
