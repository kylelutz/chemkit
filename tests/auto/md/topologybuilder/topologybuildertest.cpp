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

#include "topologybuildertest.h"

#include <chemkit/atom.h>
#include <chemkit/molecule.h>
#include <chemkit/topology.h>
#include <chemkit/topologybuilder.h>

void TopologyBuilderTest::phenol()
{
    chemkit::Molecule phenol("c1ccccc1O", "smiles");
    QCOMPARE(phenol.formula(), std::string("C6H6O"));

    chemkit::TopologyBuilder builder;
    builder.setAtomTyper("uff");
    builder.addMolecule(&phenol);

    boost::shared_ptr<chemkit::Topology> topology = builder.topology();
    QVERIFY(topology != 0);

    QCOMPARE(topology->size(), size_t(13));
    QCOMPARE(topology->bondedInteractionCount(), size_t(13));

    for(size_t i = 0; i < phenol.size(); i++){
        const chemkit::Atom *atom = phenol.atom(i);

        if(atom->is(chemkit::Atom::Hydrogen)){
            QCOMPARE(topology->type(i), std::string("H_"));
        }
        else if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(topology->type(i), std::string("C_R"));
        }
        else if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(topology->type(i), std::string("O_3"));
        }
    }
}

QTEST_APPLESS_MAIN(TopologyBuilderTest)
