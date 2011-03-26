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

#include "inchilineformat.h"

#include "../../3rdparty/inchi/inchi_api.h"

InchiLineFormat::InchiLineFormat()
    : chemkit::LineFormat("inchi")
{
}

bool InchiLineFormat::read(const std::string &formula, chemkit::Molecule *molecule)
{
    // verify formula
    if(formula.empty()){
        return 0;
    }

    QString formulaString = formula.c_str();

    // add `InChI=` to the start if it is not there
    if(formula.compare(0, 6, "InChI=") != 0){
        formulaString.prepend("InChI=");
    }

    // setup input struct
    inchi_InputINCHI input;

    input.szInChI = new char[formulaString.size()+1];
    strncpy(input.szInChI, formulaString.toAscii().data(), formulaString.size()+1);
    input.szOptions = 0;

    // get inchi output
    inchi_OutputStruct output;
    int ret = GetStructFromStdINCHI(&input, &output);
    Q_UNUSED(ret);

    // build molecule from inchi output
    QVector<chemkit::Atom *> atoms(output.num_atoms);

    bool addHydrogens = option("add-hydrogens").toBool();

    // add atoms
    for(int i = 0; i < output.num_atoms; i++){
        inchi_Atom *inchiAtom = &output.atom[i];

        chemkit::Atom *atom = molecule->addAtom(inchiAtom->elname);
        atoms[i] = atom;

        // add implicit hydrogens (if enabled)
        if(addHydrogens){
            for(int j = 0; j < inchiAtom->num_iso_H[0]; j++){
                chemkit::Atom *hydrogen = molecule->addAtom(chemkit::Atom::Hydrogen);
                molecule->addBond(atom, hydrogen);
            }
        }
    }

    // add bonds
    for(int i = 0; i < output.num_atoms; i++){
        inchi_Atom *inchiAtom = &output.atom[i];

        chemkit::Atom *atom = atoms[i];

        for(int j = 0; j < inchiAtom->num_bonds; j++){
            chemkit::Atom *neighbor = atoms[inchiAtom->neighbor[j]];
            molecule->addBond(atom, neighbor, inchiAtom->bond_type[j]);
        }
    }

    // free input and output structures
    delete [] input.szInChI;
    FreeStructFromStdINCHI(&output);

    return true;
}

std::string InchiLineFormat::write(const chemkit::Molecule *molecule)
{
    // check for valid molecule
    if(molecule->atomCount() > 1024){
        setErrorString("InChI does not support molecules with more that 1024 atoms.");
        return std::string();
    }

    // setup inchi input structure
    inchi_Input input;

    input.atom = new inchi_Atom[molecule->atomCount()];
    input.stereo0D = 0;
    input.szOptions = 0;
    input.num_atoms = molecule->atomCount();
    input.num_stereo0D = 0;

    inchi_Atom *inputAtom = &input.atom[0];

    QList<const chemkit::Atom *> chiralAtoms;

    foreach(const chemkit::Atom *atom, molecule->atoms()){

        // coordinates
        inputAtom->x = 0;
        inputAtom->y = 0;
        inputAtom->z = 0;

        // bonds and neighbors
        int neighborCount = 0;
        foreach(const chemkit::Bond *bond, atom->bonds()){
            const chemkit::Atom *neighbor = bond->otherAtom(atom);

            if(neighbor->index() < atom->index())
                continue;

            inputAtom->neighbor[neighborCount] = neighbor->index();
            inputAtom->bond_type[neighborCount] = bond->order();

            inputAtom->bond_stereo[neighborCount] = INCHI_BOND_STEREO_NONE;

            neighborCount++;
        }
        inputAtom->num_bonds = neighborCount;

        // element symbol
        strncpy(inputAtom->elname, atom->symbol().c_str(), ATOM_EL_LEN);

        // isotopic hydrogens
        inputAtom->num_iso_H[0] = -1;
        inputAtom->num_iso_H[1] = 0;
        inputAtom->num_iso_H[2] = 0;
        inputAtom->num_iso_H[3] = 0;

        // misc data
        inputAtom->isotopic_mass = 0;
        inputAtom->radical = 0;
        inputAtom->charge = 0;

        // chiral atoms
        if(atom->isChiral()){
            chiralAtoms.append(atom);
        }

        // move pointer to the next position in the atom array
        inputAtom++;
    }

    // add stereochemistry if enabled
    bool stereochemistry = option("stereochemistry").toBool();

    if(stereochemistry){
        input.stereo0D = new inchi_Stereo0D[chiralAtoms.size()];
        input.num_stereo0D = chiralAtoms.size();

        int chiralIndex = 0;

        foreach(const chemkit::Atom *atom, chiralAtoms){
            inchi_Stereo0D *stereo = &input.stereo0D[chiralIndex];
            qMemSet(stereo, 0, sizeof(*stereo));
            stereo->central_atom = atom->index();

            int neighborIndex = 0;
            foreach(const chemkit::Atom *neighbor, atom->neighbors()){
                stereo->neighbor[neighborIndex++] = neighbor->index();
            }

            stereo->type = INCHI_StereoType_Tetrahedral;

            if(atom->chirality() == chemkit::Atom::R){
                stereo->parity = INCHI_PARITY_EVEN;
            }
            else if(atom->chirality() == chemkit::Atom::S){
                stereo->parity = INCHI_PARITY_ODD;
            }
            else if(atom->chirality() == chemkit::Atom::UnspecifiedChirality){
                stereo->parity = INCHI_PARITY_UNDEFINED;
            }
            else{
                stereo->parity = INCHI_PARITY_UNKNOWN;
            }

            chiralIndex++;
        }
    }
    else{
        input.stereo0D = 0;
        input.num_stereo0D = 0;
    }

    // create inchi generator object
    INCHIGEN_HANDLE generator = STDINCHIGEN_Create();

    // setup generator object
    INCHIGEN_DATA generatorData;
    int ret = STDINCHIGEN_Setup(generator, &generatorData, &input);

    // perform structure normalization
    ret = STDINCHIGEN_DoNormalization(generator, &generatorData);

    // perform structure canonicalization
    ret = STDINCHIGEN_DoCanonicalization(generator, &generatorData);

    // write inchi output structure
    inchi_Output output;
    ret = STDINCHIGEN_DoSerialization(generator, &generatorData, &output);

    // get inchi string from output
    std::string inchiString;
    if(output.szInChI){
        inchiString = output.szInChI;
    }

    // destroy inchi input structure
    delete [] input.atom;
    delete [] input.stereo0D;

    // destroy inchi generator object
    STDINCHIGEN_Destroy(generator);

    return inchiString;
}

QVariant InchiLineFormat::defaultOption(const std::string &name) const
{
    if(name == "stereochemistry")
        return QVariant(true);
    else if(name == "add-hydrogens")
        return QVariant(true);
    else
        return QVariant();
}
