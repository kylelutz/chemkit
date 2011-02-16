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

#include "protein.h"

#include "molecule.h"
#include "proteinchain.h"

namespace chemkit {

// === ProteinPrivate ====================================================== //
class ProteinPrivate
{
    public:
        Molecule *molecule;
        QList<ProteinChain *> chains;
};

// === Protein ============================================================= //
/// \class Protein protein.h chemkit/protein.h
/// \ingroup chemkit
/// \brief The Protein class represents a protein biomolecule.
///
/// \sa AminoAcid, ProteinChain, BiochemicalFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein.
Protein::Protein()
    : d(new ProteinPrivate)
{
    d->molecule = new Molecule;
}

/// Destroys a protein.
Protein::~Protein()
{
    delete d->molecule;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of amino acid residues in the protein.
int Protein::size() const
{
    return residueCount();
}

/// Returns the molecule for the protein.
Molecule* Protein::molecule() const
{
    return d->molecule;
}

// --- Structure ----------------------------------------------------------- //
/// Creates, adds, and returns a new chain to the protein.
ProteinChain* Protein::addChain()
{
    ProteinChain *chain = new ProteinChain(this);
    d->chains.append(chain);
    return chain;
}

/// Removes a chain from the protein.
void Protein::removeChain(ProteinChain *chain)
{
    d->chains.removeOne(chain);
}

/// Returns a list of all the chains in the protein.
QList<ProteinChain *> Protein::chains() const
{
    return d->chains;
}

/// Returns the protein chain at index.
ProteinChain* Protein::chain(int index) const
{
    return d->chains.value(index, 0);
}

/// Returns the number of chains in the protein.
int Protein::chainCount() const
{
    return chains().size();
}

/// Returns a list of all the amino acid residues in the protein.
QList<AminoAcid *> Protein::residues() const
{
    QList<AminoAcid *> residues;

    foreach(const ProteinChain *chain, d->chains){
        residues.append(chain->residues());
    }

    return residues;
}

/// Returns the number of residues in the protein.
int Protein::residueCount() const
{
    int count = 0;

    foreach(const ProteinChain *chain, d->chains){
        count += chain->residueCount();
    }

    return count;
}

} // end chemkit namespace
