/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "similaritydescriptor.h"

#include "fingerprint.h"

namespace chemkit {

// === SimilarityDescriptorPrivate ========================================= //
class SimilarityDescriptorPrivate
{
public:
    const Molecule *molecule;
    Fingerprint *fingerprint;
};

// === SimilarityDescriptor ================================================ //
/// \class SimilarityDescriptor similaritydescriptor.h chemkit/similaritydescriptor.h
/// \ingroup chemkit
/// \brief The SimilarityDescriptor class is a similarity molecular
///        descriptor.
///
/// The similarity descriptor class is a molecular descriptor which
/// returns a value representing the similarity between molecules.
/// By default, the similarity value is the tanimoto coefficent
/// between the FP2 fingerprints for each molecule.
///
/// \see MolecularDescriptor, Fingerprint

// --- Construction and Destruction ---------------------------------------- //
/// Create a new similarity descriptor.
SimilarityDescriptor::SimilarityDescriptor(const Molecule *molecule,
                                           const std::string &fingerprint)
    : MolecularDescriptor("similarity"),
      d(new SimilarityDescriptorPrivate)
{
    d->molecule = molecule;
    d->fingerprint = Fingerprint::create(fingerprint);
}

/// Destroys the similarity descriptor object.
SimilarityDescriptor::~SimilarityDescriptor()
{
    delete d->fingerprint;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule to \p molecule.
void SimilarityDescriptor::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the molecule.
const Molecule* SimilarityDescriptor::molecule() const
{
    return d->molecule;
}

/// Sets the fingerprint to \p name.
void SimilarityDescriptor::setFingerprint(const std::string &name)
{
    delete d->fingerprint;
    d->fingerprint = Fingerprint::create(name);
}

/// Returns the fingerprint used to calculate similarity.
std::string SimilarityDescriptor::fingerprint() const
{
    if(d->fingerprint){
        return d->fingerprint->name();
    }

    return std::string();
}

// --- Descriptor ---------------------------------------------------------- //
/// Returns the similarity to \p molecule.
Variant SimilarityDescriptor::value(const Molecule *molecule) const
{
    if(!d->molecule || !d->fingerprint){
        return 0;
    }

    Bitset a = d->fingerprint->value(d->molecule);
    Bitset b = d->fingerprint->value(molecule);

    return Fingerprint::tanimotoCoefficient(a, b);
}

} // end chemkit namespace
