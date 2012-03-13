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
size_t Residue::size() const
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
    }
}

/// Returns a list of all the atoms in the residue.
std::vector<Atom *> Residue::atoms() const
{
    return d->atoms;
}

/// Returns the number of atoms in the residue.
size_t Residue::atomCount() const
{
    return d->atoms.size();
}

/// Returns a list of all the bonds in the residue.
std::vector<Bond *> Residue::bonds() const
{
    std::vector<Bond *> bonds;

    for(size_t i = 0; i < atomCount(); i++){
        for(size_t j = i + 1; j < atomCount(); j++){
            Bond *bond = d->atoms[i]->bondTo(d->atoms[j]);

            if(bond){
                bonds.push_back(bond);
            }
        }
    }

    return bonds;
}

/// Returns the number of bonds in the residue.
size_t Residue::bondCount() const
{
    return bonds().size();
}

/// Returns \c true if the residue contains the atom.
bool Residue::contains(const Atom *atom) const
{
    if(std::find(d->atoms.begin(), d->atoms.end(), atom) != d->atoms.end()){
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
