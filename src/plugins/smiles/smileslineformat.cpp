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

// This file contains the SmilesLineFormat class which implements
// reading and writing of SMILES strings.
//
// References:
//   OpenSMILES: http://www.opensmiles.org
//   Daylight Theory Manual: http://www.daylight.com/dayhtml/doc/theory/index.html
//   Original Paper: [Weininger 1988]

#include "smileslineformat.h"

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
    else{
        int atomicNumber = chemkit::Element::atomicNumber(toupper(first));
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
    Q_UNUSED(p);

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
QVariant SmilesLineFormat::defaultOption(const QString &name) const
{
    if(name == "stereochemistry")
        return QVariant(true);
    else if(name == "add-hydrogens")
        return QVariant(true);
    else if(name == "kekulize")
        return QVariant(false);
    else
        return QVariant();
}

// --- Input and Output ---------------------------------------------------- //
bool SmilesLineFormat::read(const QString &formula, chemkit::Molecule *molecule)
{
    return read(formula.toAscii().constData(), molecule);
}

bool SmilesLineFormat::read(const char *formula, chemkit::Molecule *molecule)
{
    const char *p = formula;
    int number = 0;
    chemkit::Atom *atom = 0;
    chemkit::Bond *bond = 0;
    chemkit::Atom *lastAtom = 0;
    int bondOrder = 1;
    QHash<chemkit::Atom *, int> charges; // atom -> formal charge
    QList<chemkit::Atom *> organicAtoms;
    QList<chemkit::Bond *> aromaticBonds;
    bool aromatic = false;
    BranchState branchState;
    QStack<BranchState> branchRoots;
    RingState ringState;
    QHash<int, RingState> rings;

    // go to initial state
    if(isTerminator(*p)) goto done;
    else if(islower(*p)) goto aromatic_atom;
    else if(*p == '[')   goto bracket_atom;
    else if(isupper(*p)) goto organic_atom;
    else                 goto parse_error;

bracket_atom:
    Q_ASSERT(*p == '[');
    p++; // move past opening bracket

    // mass number
    if(isdigit(*p)){
        number = readNumber(&p);
    }

    // symbol
    if(isupper(*p)){
        if(islower(*(p+1))){
            atom = molecule->addAtom(chemkit::Element::atomicNumber(p, 2));
            p += 2;
        }
        else{
            atom = molecule->addAtom(chemkit::Element::atomicNumber(*p));
            p++; // move past atom symbol
        }

        aromatic = false;
    }
    else if(islower(*p)){
        atom = molecule->addAtom(readAromaticSymbol(&p));
        aromatic = true;
        organicAtoms.append(atom);
    }
    else{
        goto parse_error;
    }

    if(!atom){
        goto invalid_atom_error;
    }

    // add bond to last atom
    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);
        }
        bondOrder = chemkit::Bond::Single;

        if(aromatic){
            aromaticBonds.append(bond);
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
            atom->setChirality(chemkit::Atom::S);
            p++; // move past second chirality symbol
        }
        else{
            atom->setChirality(chemkit::Atom::R);
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
    if(!atom)
        goto invalid_atom_error;

    organicAtoms.append(atom);

    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);
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
    Q_ASSERT(islower(*p));

    atom = molecule->addAtom(readAromaticSymbol(&p));
    if(!atom){
        goto invalid_atom_error;
    }

    organicAtoms.append(atom);

    if(lastAtom){
        if(bondOrder){
            bond = molecule->addBond(atom, lastAtom, bondOrder);

            if(aromatic){
                aromaticBonds.append(bond);
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

    if(rings.contains(number)){
        // ring closure
        ringState = rings.take(number);
        bond = molecule->addBond(ringState.firstAtom, lastAtom, ringState.bondOrder);

        if(aromatic && ringState.aromatic){
            aromaticBonds.append(bond);
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
    Q_ASSERT(*p == '(');
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
    Q_ASSERT(*p == ')');

    if(branchRoots.isEmpty())
        goto parse_error;

    p++; // move past closing parenthesis

    // restore state to what it was before the branch
    branchState = branchRoots.pop();
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
    setErrorString(QString("Error parsing smiles at character #%1 ('%2').")
                       .arg(p - formula)
                       .arg(*p));
    goto error;

invalid_atom_error:
    setErrorString(QString("Invalid atom symbol at character #%1 ('%2').")
                        .arg(p - formula)
                        .arg(*p));
    goto error;

error:
    return false;

done:
    // kekulize aromatic bonds
    if(!aromaticBonds.isEmpty()){
        Kekulizer::kekulize(aromaticBonds);
    }

    // add hydrogens (if enabled)
    if(option("add-hydrogens").toBool()){
        foreach(chemkit::Atom *atom, organicAtoms){
            while(atom->formalCharge() < 0){
                chemkit::Atom *hydrogen = molecule->addAtom(chemkit::Atom::Hydrogen);
                molecule->addBond(atom, hydrogen);
            }
        }
    }

    return true;
}

QString SmilesLineFormat::write(const chemkit::Molecule *molecule)
{
    bool kekulize = option("kekulize").toBool();

    return SmilesGraph(molecule).toString(kekulize);
}
