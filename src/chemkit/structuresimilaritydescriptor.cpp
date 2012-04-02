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

#include "structuresimilaritydescriptor.h"

#include <boost/make_shared.hpp>

#include "atom.h"
#include "molecule.h"
#include "substructurequery.h"

namespace chemkit {

// === StructureSimilarityDescriptorPrivate ================================ //
class StructureSimilarityDescriptorPrivate
{
public:
    boost::shared_ptr<Molecule> molecule;
};

// === StructureSimilarityDescriptor ======================================= //
/// \class StructureSimilarityDescriptor structuresimilaritydescriptor.h chemkit/structuresimilaritydescriptor.h
/// \ingroup chemkit
/// \brief The StructureSimilarityDescriptor class calculates similarity
///        between molecules based on their structures.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new structure similarity descriptor.
StructureSimilarityDescriptor::StructureSimilarityDescriptor()
    : MolecularDescriptor("structure-similarity"),
      d(new StructureSimilarityDescriptorPrivate)
{
}

/// Creates a new structure similarity descriptor with \p molecule.
StructureSimilarityDescriptor::StructureSimilarityDescriptor(const boost::shared_ptr<Molecule> &molecule)
    : MolecularDescriptor("structure-similarity"),
      d(new StructureSimilarityDescriptorPrivate)
{
    d->molecule = molecule;
}

/// Destroys the structure similarity descriptor object.
StructureSimilarityDescriptor::~StructureSimilarityDescriptor()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the query molecule for the descriptor to \p molecule.
void StructureSimilarityDescriptor::setMolecule(const boost::shared_ptr<Molecule> &molecule)
{
    d->molecule = molecule;
}

/// Returns the query molecule for the descriptor.
boost::shared_ptr<Molecule> StructureSimilarityDescriptor::molecule() const
{
    return d->molecule;
}

// --- Descriptor ---------------------------------------------------------- //
/// Returns the structure similarity value for \p molecule.
Variant StructureSimilarityDescriptor::value(const Molecule *molecule) const
{
    if(!d->molecule){
        return 0;
    }

    SubstructureQuery query(d->molecule);

    size_t a = d->molecule->atomCount() - d->molecule->atomCount(Atom::Hydrogen);
    size_t b = molecule->atomCount() - molecule->atomCount(Atom::Hydrogen);
    size_t c = query.maximumMapping(molecule).size();

    return Real(c) / Real(a + b - c);
}

} // end chemkit namespace
