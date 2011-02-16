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

#include "proteinchain.h"

#include "protein.h"

namespace chemkit {

// === ProteinChainPrivate ================================================= //
class ProteinChainPrivate
{
    public:
        Protein *protein;
        QList<AminoAcid *> residues;
};

// === ProteinChain ======================================================== //
/// \class ProteinChain proteinchain.h chemkit/proteinchain.h
/// \ingroup chemkit
/// \brief The ProteinChain class represents a single chain of amino
///        acid residues in a protein.
///
/// ProteinChain objects are created with the Protein::addChain()
/// method and destroyed with the Protein::removeChain() method.
///
/// \see Protein, AminoAcid

// --- Construction and Destruction ---------------------------------------- //
ProteinChain::ProteinChain(Protein *protein)
    : d(new ProteinChainPrivate)
{
    d->protein = protein;
}

ProteinChain::~ProteinChain()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of amino acid residues in the chain.
int ProteinChain::size() const
{
    return residueCount();
}

/// Returns the protein the chain is a part of.
Protein* ProteinChain::protein() const
{
    return d->protein;
}

/// Returns the molecule for the protein.
Molecule* ProteinChain::molecule() const
{
    return d->protein->molecule();
}

// --- Structure ----------------------------------------------------------- //
/// Adds an amino acid residue to the end of the chain.
void ProteinChain::addResidue(AminoAcid *residue)
{
    d->residues.append(residue);
}

/// Creates a new amino acid residue and adds it to the end of the
/// chain.
AminoAcid* ProteinChain::addNewResidue()
{
    AminoAcid* residue = new AminoAcid(molecule());
    molecule()->addResidue(residue);

    addResidue(residue);

    return residue;
}

/// Removes the residue from the chain.
void ProteinChain::removeResidue(AminoAcid *residue)
{
    d->residues.removeOne(residue);
}

/// Returns a list of all the amino acid residues in the chain.
QList<AminoAcid *> ProteinChain::residues() const
{
    return d->residues;
}

/// Returns the amino acid residue at index.
AminoAcid* ProteinChain::residue(int index) const
{
    return d->residues.value(index, 0);
}

/// Returns the number of amino acid residues in the chain.
int ProteinChain::residueCount() const
{
    return d->residues.size();
}

/// Returns the number of amino acid residues of the given type in the
/// chain.
int ProteinChain::residueCount(AminoAcid::AminoAcidType type) const
{
    int count = 0;

    foreach(AminoAcid *aminoAcid, d->residues)
        if(aminoAcid->type() == type)
            count++;

    return count;
}

/// Returns the amino acid sequence as a string of one letter symbols.
QString ProteinChain::sequenceString() const
{
    QString string;

    foreach(AminoAcid *residue, d->residues){
        string.append(residue->letter());
    }

    return string;
}

/// Returns the sequence number of the residue.
int ProteinChain::sequenceNumber(const AminoAcid *residue) const
{
    return d->residues.indexOf(const_cast<AminoAcid *>(residue)) + 1;
}

} // end chemkit namespace
