/******************************************************************************
**
** Copyright (C) 2012 Kitware, Inc.
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

#include "gravitationalindexdescriptor.h"

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/geometry.h>
#include <chemkit/molecule.h>

// === GravitationalIndexDescriptor ======================================== //
GravitationalIndexDescriptor::GravitationalIndexDescriptor()
    : chemkit::MolecularDescriptor("gravitational-index")
{
    setDimensionality(3);
}

chemkit::Variant GravitationalIndexDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real value = 0;

    for(size_t i = 0; i < molecule->atomCount(); i++){
        const chemkit::Atom *a = molecule->atom(i);

        for(size_t j = i + 1; j < molecule->atomCount(); j++){
            const chemkit::Atom *b = molecule->atom(j);

            chemkit::Real r2 = chemkit::geometry::distanceSquared(a->position(),
                                                                  b->position());

            value += (a->mass() * b->mass()) / r2;
        }
    }

    return value;
}

// === BondedGravitationalIndexDescriptor ================================== //
BondedGravitationalIndexDescriptor::BondedGravitationalIndexDescriptor()
    : chemkit::MolecularDescriptor("bonded-gravitational-index")
{
    setDimensionality(3);
}

chemkit::Variant BondedGravitationalIndexDescriptor::value(const chemkit::Molecule *molecule) const
{
    chemkit::Real value = 0;

    foreach(const chemkit::Bond *bond, molecule->bonds()){
        const chemkit::Atom *a = bond->atom1();
        const chemkit::Atom *b = bond->atom2();

        chemkit::Real r2 = chemkit::geometry::distanceSquared(a->position(),
                                                              b->position());

        value += (a->mass() * b->mass()) / r2;
    }

    return value;
}
