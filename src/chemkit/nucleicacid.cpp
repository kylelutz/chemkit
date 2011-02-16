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

#include "nucleicacid.h"

#include "molecule.h"
#include "nucleicacidchain.h"

namespace chemkit {

// === NucleicAcidPrivate ================================================== //
class NucleicAcidPrivate
{
    public:
        Molecule *molecule;
        QList<NucleicAcidChain *> chains;
};

// === NucleicAcid ========================================================= //
/// \class NucleicAcid nucleicacid.h chemkit/nucleicacid.h
/// \ingroup chemkit
/// \brief The NucleicAcid class represents a nucleic acid
///        biomolecule (DNA or RNA).

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new nucleic acid.
NucleicAcid::NucleicAcid()
    : d(new NucleicAcidPrivate)
{
    d->molecule = new Molecule;
}

/// Destroys the nucleic acid.
NucleicAcid::~NucleicAcid()
{
    delete d->molecule;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of residues in the nucleic acid.
int NucleicAcid::size() const
{
    return residueCount();
}

/// Returns the molecule for the nucleic acid.
Molecule* NucleicAcid::molecule() const
{
    return d->molecule;
}

// --- Structure ----------------------------------------------------------- //
/// Adds a new chain to the nucleic acid.
NucleicAcidChain* NucleicAcid::addChain()
{
    NucleicAcidChain *chain = new NucleicAcidChain(this);
    d->chains.append(chain);
    return chain;
}

/// Removes a chain from the nucleic acid.
void NucleicAcid::removeChain(NucleicAcidChain *chain)
{
    d->chains.removeOne(chain);
}

/// Returns the nucleic acid chain at \p index.
NucleicAcidChain* NucleicAcid::chain(int index) const
{
    return d->chains.value(index, 0);
}

/// Returns a list of all the chains in the nucleic acid.
QList<NucleicAcidChain *> NucleicAcid::chains() const
{
    return d->chains;
}

/// Returns a list of all the chains in the nucleic acid.
int NucleicAcid::chainCount() const
{
    return d->chains.size();
}

/// Returns the number of residues in the nucleic acid.
int NucleicAcid::residueCount() const
{
    return 0;
}

} // end chemkit namespace
