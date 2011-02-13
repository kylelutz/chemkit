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

#include <QtCore>

#include <chemkit/molecule.h>
#include <chemkit/forcefield.h>
#include <chemkit/chemicalfile.h>

int main(int argc, char *argv[])
{
    QTextStream out(stdout);
    QTextStream err(stderr);

    if(argc < 2){
        err << "Usage: " << argv[0] << " FILENAME" << "\n";
        return -1;
    }

    QString fileName = argv[1];

    chemkit::ChemicalFile file(fileName);
    bool ok = file.read();
    if(!ok){
        err << "Failed to read file: " << fileName << "\n";
        return -1;
    }

    chemkit::Molecule *molecule = file.molecule();
    if(!molecule){
        err << "File contains no molecules.\n";
        return -1;
    }

    chemkit::ForceField *uff = chemkit::ForceField::create("uff");
    if(!uff){
        err << "UFF force field plugin not found.\n";
        return -1;
    }

    uff->addMolecule(molecule);
    uff->setup();

    if(!uff->isSetup()){
        err << "Failed to parameterize force field.\n";
        delete uff;
        return -1;
    }

    chemkit::Float energy = uff->energy();

    out << "Formula: " << molecule->formula() << "\n";
    out << "Energy: " << energy << " kcal/mol\n";

    delete uff;

	return 0;
}
