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

#include "residue.h"

#include <map>
#include <algorithm>

#include "atom.h"
#include "molecule.h"

namespace chemkit {

// === ResiduePrivate ====================================================== //
class ResiduePrivate
{
    public:
        int type;
        Molecule *molecule;
        std::vector<Atom *> atoms;
        std::map<std::string, const Atom *> types;
};

// === Residue ============================================================= //
/// \class Residue residue.h chemkit/residue.h
/// \ingroup chemkit
/// \brief The Residue class represents a single monomer in a larger
///        molecule.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new residue of given type.
Residue::Residue(Molecule *molecule, int type)
    : d(new ResiduePrivate)
{
    d->type = type;
    d->molecule = molecule;
}

/// Destroys the residue.
Residue::~Residue()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the residue type.
int Residue::residueType() const
{
    return d->type;
}

/// Returns the number of atoms in the residue.
int Residue::size() const
{
    return atomCount();
}

/// Returns a one letter symbol for the residue;
char Residue::letter() const
{
    return ' ';
}

/// Returns the molecule the residue is a part of.
Molecule* Residue::molecule() const
{
    return d->molecule;
}

// --- Structure ----------------------------------------------------------- //
/// Adds an atom to the residue.
void Residue::addAtom(Atom *atom)
{
    if(atom->molecule() != d->molecule || contains(atom)){
        return;
    }

    d->atoms.push_back(atom);
}

/// Removes an atom from the residue.
void Residue::removeAtom(Atom *atom)
{
    std::vector<Atom *>::iterator location = std::find(d->atoms.begin(), d->atoms.end(), atom);
    if(location != d->atoms.end()){
        d->atoms.erase(location);

        if(atom->residue() == this){
            atom->setResidue(0);
        }
    }
}

/// Returns a list of all the atoms in the residue.
std::vector<Atom *> Residue::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the residue.
int Residue::atomCount() const
{
    return d->atoms.size();
}

/// Returns a list of all the bonds in the residue.
std::vector<Bond *> Residue::bonds() const
{
    std::vector<Bond *> bonds;

    for(int i = 0; i < atomCount(); i++){
        for(int j = i + 1; j < atomCount(); j++){
            Bond *bond = d->atoms[i]->bondTo(d->atoms[j]);

            if(bond){
                bonds.push_back(bond);
            }
        }
    }

    return bonds;
}

/// Returns the number of bonds in the residue.
int Residue::bondCount() const
{
    return bonds().size();
}

/// Returns \c true if the residue contains the atom.
bool Residue::contains(const Atom *atom) const
{
    if(atom->residue() == this){
        return true;
    }
    else if(std::find(d->atoms.begin(), d->atoms.end(), atom) != d->atoms.end()){
        return true;
    }
    else{
        return false;
    }
}

/// Returns \c true if the residue contains the bond.
bool Residue::contains(const Bond *bond) const
{
    const std::vector<Bond *> &bonds = this->bonds();

    return std::find(bonds.begin(), bonds.end(), bond) != bonds.end();
}

// --- Atom Types ---------------------------------------------------------- //
/// Sets the type for atom in the residue.
void Residue::setAtomType(const Atom *atom, const std::string &type)
{
    d->types[type] = atom;
}

/// Returns the type for atom in the residue.
std::string Residue::atomType(const Atom *atom) const
{
    for(std::map<std::string, const Atom *>::iterator i = d->types.begin(); i != d->types.end(); ++i){
        if(i->second == atom){
            return i->first;
        }
    }

    return std::string();
}

/// Returns the atom with type or 0 if no atom has type.
Atom* Residue::atom(const std::string &type) const
{
    return const_cast<Atom *>(d->types[type]);
}

} // end chemkit namespace
