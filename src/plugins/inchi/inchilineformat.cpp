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

#include "inchilineformat.h"

#include <QString>
#include <QVector>

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

chemkit::Variant InchiLineFormat::defaultOption(const std::string &name) const
{
    if(name == "stereochemistry")
        return true;
    else if(name == "add-hydrogens")
        return true;
    else
        return chemkit::Variant();
}
