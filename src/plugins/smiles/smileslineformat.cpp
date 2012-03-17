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

// This file contains the SmilesLineFormat class which implements
// reading and writing of SMILES strings.
//
// References:
//   OpenSMILES: http://www.opensmiles.org
//   Daylight Theory Manual: http://www.daylight.com/dayhtml/doc/theory/index.html
//   Original Paper: [Weininger 1988]

#include "smileslineformat.h"

#include <stack>
#include <vector>

#include <boost/format.hpp>

#include <chemkit/foreach.h>

#include "smiles.h"
#include "kekulizer.h"
#include "smilesgraph.h"

namespace {

// Returns true if the character represents a bond symbol.
bool isBond(char c)
{
    return c == '-' || c == '=' || c == '#' || c == '$' ||
           c == '.' || c == '/' || c == '\\';
}

// Returns true if the character represents the end of a SMILES
// string (e.g. the null terminator or a whitespace character).
bool isTerminator(char c)
{
    return c == '\0' || isspace(c);
}

// Returns true if the character represents a ring identifier.
bool isRing(char c)
{
    return isdigit(c) || c == '%';
}

int readAromaticSymbol(const char **p)
{
    char first = **p;
    char second = *(*p+1);

    (*p)++;

    if(first == 'a' && second == 's'){
        (*p)++;
        return chemkit::Atom::Arsenic;
    }
    else if(first == 's' && second == 'e'){
        (*p)++;
        return chemkit::Atom::Selenium;
    }
    else if(first == 't' && second == 'e'){
        (*p)++;
        return chemkit::Atom::Tellurium;
    }
    else{
        int atomicNumber = chemkit::Element::fromSymbol(toupper(first)).atomicNumber();
        if(!atomicNumber){
            (*p)--;
            return 0;
        }
        else{
            return atomicNumber;
        }
    }
}

int readOrganicSymbol(const char **p)
{
    char first = **p;
    char second = *(*p+1);

    (*p)++;

    switch(first){
        case 'B':
            if(second == 'r'){
                (*p)++;
                return chemkit::Atom::Bromine;
            }
            else{
                return chemkit::Atom::Boron;
            }
        case 'C':
            if(second == 'l'){
                (*p)++;
                return chemkit::Atom::Chlorine;
            }
            else{
                return chemkit::Atom::Carbon;
            }
        case 'N':
            return chemkit::Atom::Nitrogen;
        case 'O':
            return chemkit::Atom::Oxygen;
        case 'P':
            return chemkit::Atom::Phosphorus;
        case 'S':
            return chemkit::Atom::Sulfur;
        case 'F':
            return chemkit::Atom::Fluorine;
        case 'I':
            return chemkit::Atom::Iodine;
        default:
            (*p)--;
            return 0;
    }
}

int readNumber(const char **p)
{
    int number = *(*p)++ - '0';

    while(isdigit(**p)){
        number = number * 10 + *(*p)++ - '0';
    }

    return number;
}

int readCharge(const char **p)
{
    CHEMKIT_UNUSED(p);

    return 0;
}

struct BranchState
{
    chemkit::Atom *lastAtom;
    int bondOrder;
    bool aromatic;
};

struct RingState
{
    chemkit::Atom *firstAtom;
    int bondOrder;
    bool aromatic;
};

} // end anonymous namespace

// === SmilesLineFormat ==================================================== //
SmilesLineFormat::SmilesLineFormat()
    : chemkit::LineFormat("smiles")
{
}

// --- Options ------------------------------------------------------------- //
chemkit::Variant SmilesLineFormat::defaultOption(const std::string &name) const
{
    if(name == "stereochemistry")
        return true;
    else if(name == "add-implicit-hydrogens")
        return true;
    else if(name == "kekulize")
        return false;
    else
        return chemkit::Variant();
}

// --- Input and Output ---------------------------------------------------- //
chemkit::Molecule* SmilesLineFormat::read(const std::string &formula)
{
    return read(formula.c_str());
}

chemkit::Molecule* SmilesLineFormat::read(const char *formula)
{
    const char *p = formula;
    int number = 0;
    chemkit::Atom *atom = 0;
    chemkit::Bond *bond = 0;
    chemkit::Atom *lastAtom = 0;
    chemkit::Bond *lastDoubleBond = 0;
    int bondOrder = 1;
    std::map<chemkit::Atom *, int> charges; // atom -> formal charge
    std::vector<chemkit::Atom *> organicAtoms;
    std::vector<chemkit::Bond *> aromaticBonds;
    bool aromatic = false;
    BranchState branchState;
    std::stack<BranchState> branchRoots;
    RingState ringState;
    std::map<int, RingState> rings;

    enum BondStereo {
        Up = 1,
        Down = 2
    };

    int bondStereo = 0;

    // create molecule
    chemkit::Molecule *molecule = new chemkit::Molecule;

    // go to initial state
    if(isTerminator(*p)) goto done;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isupper(*p)) goto organic_atom;
    else                 goto parse_error;

bracket_atom:
    assert(*p == '[');
    p++; // move past opening bracket

    // mass number
    if(isdigit(*p)){
        number = readNumber(&p);
    }

    // symbol
    if(isupper(*p)){
        if(islower(*(p+1))){
            atom = molecule->addAtom(chemkit::Element::fromSymbol(p, 2));
            p += 2;
        }
        else{
            atom = molecule->addAtom(chemkit::Element::fromSymbol(*p));
            p++; // move past atom symbol
        }

        aromatic = false;
    }
    else if(islower(*p)){
        atom = molecule->addAtom(readAromaticSymbol(&p));
        aromatic = true;
        organicAtoms.push_back(atom);
    }
    else{
        goto parse_error;
    }

    if(!atom->element().isValid()){
        goto invalid_atom_error;
    }

    // add bond to last atom
    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);
        }
        bondOrder = chemkit::Bond::Single;

        if(bond && aromatic){
            aromaticBonds.push_back(bond);
        }
    }

    lastAtom = atom;

    // set mass number
    if(number){
        atom->setMassNumber(number);
    }

    // chirality symbol
    if(*p == '@'){
        p++; // move past chirality symbol

        if(*p == '@'){
            atom->setChirality(chemkit::Stereochemistry::S);
            p++; // move past second chirality symbol
        }
        else{
            atom->setChirality(chemkit::Stereochemistry::R);
        }
    }

    // hydrogens
    if(*p == 'H'){
        p++; // move past 'H' symbol
        number = 1;

        if(isdigit(*p)){
            number = *p - '0';
            p++; // move past digit
        }

        for(int i = 0; i < number; i++){
            molecule->addBond(atom, molecule->addAtom(chemkit::Atom::Hydrogen));
        }
    }

    // charge
    if(*p == '+'){
        p++;

        if(*p == '+'){
            charges[atom] = +2;
            p++;
        }
        else if(isdigit(*p)){
            charges[atom] = readNumber(&p);
        }
        else{
            charges[atom] = +1;
        }
    }
    else if(*p == '-'){
        p++;

        if(*p == '-'){
            charges[atom] = -2;
            p++;
        }
        else if(isdigit(*p)){
            charges[atom] = -readNumber(&p);
        }
        else{
            charges[atom] = -1;
        }
    }

    // end bracket
    if(*p != ']'){
        goto parse_error;
    }

    p++; // move past closing bracket

    // go to next state
    if(isTerminator(*p)) goto done;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isBond(*p))  goto bond;
    else if(isRing(*p))  goto ring;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

organic_atom:
    atom = molecule->addAtom(readOrganicSymbol(&p));
    if(!atom->element().isValid())
        goto invalid_atom_error;

    organicAtoms.push_back(atom);

    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);

            if(bondOrder == chemkit::Bond::Double){
                lastDoubleBond = bond;
            }
        }

        bondOrder = chemkit::Bond::Single;
    }

    aromatic = false;
    lastAtom = atom;

    // go to next state
    if(isTerminator(*p)) goto done;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isBond(*p))  goto bond;
    else if(isRing(*p))  goto ring;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

aromatic_atom:
    assert(islower(*p));

    atom = molecule->addAtom(readAromaticSymbol(&p));
    if(!atom->element().isValid()){
        goto invalid_atom_error;
    }

    organicAtoms.push_back(atom);

    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);

            if(bondOrder == chemkit::Bond::Double){
                lastDoubleBond = bond;
            }

            if(aromatic){
                aromaticBonds.push_back(bond);
            }
        }
    }

    aromatic = true;
    lastAtom = atom;

    // go to next state
    if(isTerminator(*p)) goto done;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isBond(*p))  goto bond;
    else if(isRing(*p))  goto ring;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

bond:
    if(*p == '-')
        bondOrder = 1;
    else if(*p == '=')
        bondOrder = 2;
    else if(*p == '#')
        bondOrder = 3;
    else if(*p == '$')
        bondOrder = 4;
    else if(*p == '.')
        bondOrder = 0;
    else if(*p == '/')
        bondOrder = 1;
    else if(*p == '\\')
        bondOrder = 1;

    // set stereochemistry
    if(*p == '/'){
        if(bondStereo != 0 && lastDoubleBond){
            if(bondStereo == Up){
                lastDoubleBond->setStereochemistry(chemkit::Stereochemistry::E);
            }
            else if(bondStereo == Down){
                lastDoubleBond->setStereochemistry(chemkit::Stereochemistry::Z);
            }

            bondStereo = 0;
            lastDoubleBond = 0;
        }
        else{
            bondStereo = Up;
        }
    }
    else if(*p == '\\'){
        if(bondStereo != 0 && lastDoubleBond){
            if(bondStereo == Up){
                lastDoubleBond->setStereochemistry(chemkit::Stereochemistry::Z);
            }
            else if(bondStereo == Down){
                lastDoubleBond->setStereochemistry(chemkit::Stereochemistry::E);
            }

            bondStereo = 0;
            lastDoubleBond = 0;
        }
        else{
            bondStereo = Down;
        }
    }

    p++; // move past bond symbol

    // go to next state
    if(*p == '[')        goto bracket_atom;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(isRing(*p))  goto ring;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

ring:
    number = 0;

    if(*p == '%'){
        p++; // move past '%'
        number = readNumber(&p);
    }
    else{
        number = *p - '0';
        p++; // move past digit
    }

    if(rings.find(number) != rings.end()){
        // ring closure
        ringState = rings[number];
        rings.erase(number);
        bond = molecule->addBond(ringState.firstAtom, lastAtom, ringState.bondOrder);

        if(aromatic && ringState.aromatic){
            aromaticBonds.push_back(bond);
        }
    }
    else{
        // ring opening
        ringState.firstAtom = lastAtom;
        ringState.bondOrder = bondOrder;
        ringState.aromatic = aromatic;
        rings[number] = ringState;
    }

    // go to next state
    if(isTerminator(*p)) goto done;
    else if(*p == '[')   goto bracket_atom;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(isRing(*p))  goto ring;
    else if(isBond(*p))  goto bond;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

start_branch:
    assert(*p == '(');
    p++; // move past opening parenthesis

    // save current state
    branchState.lastAtom = lastAtom;
    branchState.bondOrder = bondOrder;
    branchState.aromatic = aromatic;
    branchRoots.push(branchState);

    // go to next state
    if(isupper(*p))      goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isBond(*p))  goto bond;
    else                 goto parse_error;

end_branch:
    assert(*p == ')');

    if(branchRoots.empty())
        goto parse_error;

    p++; // move past closing parenthesis

    // restore state to what it was before the branch
    branchState = branchRoots.top();
    branchRoots.pop();
    lastAtom = branchState.lastAtom;
    bondOrder = branchState.bondOrder;
    aromatic = branchState.aromatic;

    // go to next state
    if(isTerminator(*p)) goto done;
    else if(isupper(*p)) goto organic_atom;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isBond(*p))  goto bond;
    else if(isRing(*p))  goto ring;
    else if(*p == '(')   goto start_branch;
    else if(*p == ')')   goto end_branch;
    else                 goto parse_error;

parse_error:
    setErrorString((boost::format("Error parsing smiles at character #%d ('%c').") %
                       (p - formula) % *p).str());
    goto error;

invalid_atom_error:
    setErrorString((boost::format("Invalid atom symbol at character #%d ('%c').") %
                             (p - formula) % *p).str());
    goto error;

error:
    delete molecule;
    return 0;

done:
    // kekulize aromatic bonds
    if(!aromaticBonds.empty()){
        Kekulizer::kekulize(aromaticBonds);
    }

    // add implicit hydrogens (if enabled)
    if(option("add-implicit-hydrogens").toBool()){
        foreach(chemkit::Atom *atom, organicAtoms){
            while(atom->formalCharge() < 0){
                chemkit::Atom *hydrogen = molecule->addAtom(chemkit::Atom::Hydrogen);
                molecule->addBond(atom, hydrogen);
            }
        }
    }

    return molecule;
}

std::string SmilesLineFormat::write(const chemkit::Molecule *molecule)
{
    bool kekulize = option("kekulize").toBool();

    return SmilesGraph(molecule).toString(kekulize);
}
