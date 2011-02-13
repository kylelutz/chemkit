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

#include "nucleicacidchain.h"

#include "nucleotide.h"

namespace chemkit {

// === NucleicAcidChainPrivate ============================================= //
class NucleicAcidChainPrivate
{
    public:
        NucleicAcid *nucleicAcid;
        QList<Nucleotide *> residues;
};

// === NucleicAcidChain ==================================================== //
/// \class NucleicAcidChain nucleicacidchain.h chemkit/nucleicacidchain.h
/// \ingroup chemkit
/// \brief The NucleicAcidChain class represents a single chain of
///        nucleotides in a nucleic acid.
///
/// NucleicAcidChain objects are created with the
/// NucleicAcid::addChain() method and destroyed with the
/// NucleicAcid::removeChain() method.
///
/// \see NucleicAcid, Nucleotide

// --- Construction and Destruction ---------------------------------------- //
NucleicAcidChain::NucleicAcidChain(NucleicAcid *nucleicAcid)
    : d(new NucleicAcidChainPrivate)
{
    d->nucleicAcid = nucleicAcid;
}

NucleicAcidChain::~NucleicAcidChain()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the nucleic acid the chain is a part of.
NucleicAcid* NucleicAcidChain::nucleicAcid()
{
    return d->nucleicAcid;
}

/// \overload
const NucleicAcid* NucleicAcidChain::nucleicAcid() const
{
    return d->nucleicAcid;
}

// --- Structure ----------------------------------------------------------- //
/// Adds a nucleotide residue to the end of the chain.
void NucleicAcidChain::addResidue(Nucleotide *residue)
{
    d->residues.append(residue);
}

/// Removes a nucleotide residue from the chain.
void NucleicAcidChain::removeResidue(Nucleotide *residue)
{
    d->residues.removeOne(residue);
}

/// Returns a list of all the nucleotide residues in the chain.
QList<Nucleotide *> NucleicAcidChain::residues()
{
    return d->residues;
}

/// \overload
QList<const Nucleotide *> NucleicAcidChain::residues() const
{
    QList<const Nucleotide *> residues;

    foreach(const Nucleotide *residue, d->residues){
        residues.append(residue);
    }

    return residues;
}

/// Returns the number of nucleotide residues in the chain.
int NucleicAcidChain::residueCount() const
{
    return d->residues.size();
}

/// Returns \c true if the chain contains \p residue.
bool NucleicAcidChain::contains(const Nucleotide *residue) const
{
    return d->residues.contains(const_cast<Nucleotide *>(residue));
}

/// Returns a string containing the sequence of the nucleic acid
/// chain.
QString NucleicAcidChain::sequenceString() const
{
    QString string;

    foreach(const Nucleotide *residue, d->residues){
        string += residue->letter();
    }

    return string;
}

} // end chemkit namespace
