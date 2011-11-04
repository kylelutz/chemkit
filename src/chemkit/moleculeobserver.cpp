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

#include "moleculeobserver.h"

namespace chemkit {

// === MoleculeObserverPrivate ============================================= //
class MoleculeObserverPrivate
{
public:
    const Molecule *molecule;
};

// === MoleculeObserver ==================================================== //
/// \class MoleculeObserver moleculeobserver.h chemkit/moleculeobserver.h
/// \ingroup chemkit
/// \brief The MoleculeObserver class monitors a molecule for changes.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecule observer object.
MoleculeObserver::MoleculeObserver(const Molecule *molecule)
    : d(new MoleculeObserverPrivate)
{
    d->molecule = 0;

    setMolecule(molecule);
}

/// Destroys the molecule observer object.
MoleculeObserver::~MoleculeObserver()
{
    setMolecule(0);

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the observer to watch.
void MoleculeObserver::setMolecule(const Molecule *molecule)
{
    if(molecule == d->molecule){
        return;
    }

    if(d->molecule){
        d->molecule->removeObserver(this);
    }

    d->molecule = molecule;

    if(d->molecule){
        molecule->addObserver(this);
    }
}

/// Returns the molecule that the observer is watching.
const Molecule* MoleculeObserver::molecule() const
{
    return d->molecule;
}

// --- Events -------------------------------------------------------------- //
/// This virtual method is called when an atom in the molecule
/// changes.
void MoleculeObserver::atomChanged(const Atom *atom, Molecule::ChangeType changeType)
{
    CHEMKIT_UNUSED(atom);
    CHEMKIT_UNUSED(changeType);
}

/// This virtual method is called when a bond in the molecule
/// changes.
void MoleculeObserver::bondChanged(const Bond *bond, Molecule::ChangeType changeType)
{
    CHEMKIT_UNUSED(bond);
    CHEMKIT_UNUSED(changeType);
}

/// This virtual method is called when a conformer in the molecule
/// changes.
void MoleculeObserver::conformerChanged(const Conformer *conformer, Molecule::ChangeType changeType)
{
    CHEMKIT_UNUSED(conformer);
    CHEMKIT_UNUSED(changeType);
}

/// This virtual method is called when a property of the molecule
/// changes.
void MoleculeObserver::moleculeChanged(const Molecule *molecule, Molecule::ChangeType changeType)
{
    CHEMKIT_UNUSED(molecule);
    CHEMKIT_UNUSED(changeType);
}

} // end chemkit namespace

