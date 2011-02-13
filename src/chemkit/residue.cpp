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

#include "atom.h"
#include "molecule.h"

namespace chemkit {

// === ResiduePrivate ====================================================== //
class ResiduePrivate
{
    public:
        int type;
        Molecule *molecule;
        QList<Atom *> atoms;
        QHash<QString, const Atom *> types;
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

/// Returns the molecule the residue is a part of.
Molecule* Residue::molecule()
{
    return d->molecule;
}

/// \overload
const Molecule* Residue::molecule() const
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

    d->atoms.append(atom);
}

/// Removes an atom from the residue.
void Residue::removeAtom(Atom *atom)
{
    if(contains(atom)){
        d->atoms.removeOne(atom);

        if(atom->residue() == this){
            atom->setResidue(0);
        }
    }
}

/// Returns a list of all the atoms in the residue.
QList<Atom *> Residue::atoms()
{
    return d->atoms;
}

/// \overload
QList<const Atom *> Residue::atoms() const
{
    QList<const Atom *> atoms;

    foreach(const Atom *atom, const_cast<Residue *>(this)->atoms()){
        atoms.append(atom);
    }

    return atoms;
}

/// Returns the number of atoms in the residue.
int Residue::atomCount() const
{
    return d->atoms.size();
}

/// Returns a list of all the bonds in the residue.
QList<Bond *> Residue::bonds()
{
    QList<Bond *> bonds;

    for(int i = 0; i < atomCount(); i++){
        for(int j = i+1; j < atomCount(); j++){
            Bond *bond = d->atoms[i]->bondTo(d->atoms[j]);

            if(bond){
                bonds.append(bond);
            }
        }
    }

    return bonds;
}

/// \overload
QList<const Bond *> Residue::bonds() const
{
    QList<const Bond *> bonds;

    foreach(const Bond *bond, const_cast<Residue *>(this)->bonds()){
        bonds.append(bond);
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
    else if(atoms().contains(atom)){
        return true;
    }
    else{
        return false;
    }
}

/// Returns \c true if the residue contains the bond.
bool Residue::contains(const Bond *bond) const
{
    return bonds().contains(bond);
}

// --- Atom Types ---------------------------------------------------------- //
/// Sets the type for atom in the residue.
void Residue::setAtomType(const Atom *atom, const QString &type)
{
    d->types[type] = atom;
}

/// Returns the type for atom in the residue.
QString Residue::atomType(const Atom *atom) const
{
    return d->types.key(atom, QString());
}

/// Returns the atom with type or 0 if no atom has type.
Atom* Residue::atom(const QString &type)
{
    return const_cast<Atom *>(d->types.value(type, 0));
}

/// \overload
const Atom* Residue::atom(const QString &type) const
{
    return d->types.value(type, 0);
}

} // end chemkit namespace
