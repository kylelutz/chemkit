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

#include "moleculewatcher.h"

#include "molecule.h"

namespace chemkit {

// === MoleculeWatcherPrivate ============================================== //
class MoleculeWatcherPrivate
{
public:
    const Molecule *molecule;
};

// === MoleculeWatcher ===================================================== //
/// \class MoleculeWatcher moleculewatcher.h chemkit/moleculewatcher.h
/// \ingroup chemkit
/// \brief The MoleculeWatcher class monitors a molecule and emits
///        signals when changes occur.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecule watcher that monitors molecule.
MoleculeWatcher::MoleculeWatcher(const Molecule *molecule)
    : d(new MoleculeWatcherPrivate)
{
    d->molecule = 0;

    setMolecule(molecule);
}

/// Destroys the molecule watcher object.
MoleculeWatcher::~MoleculeWatcher()
{
    setMolecule(0);

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule for the watcher to monitor.
void MoleculeWatcher::setMolecule(const Molecule *molecule)
{
    if(molecule == d->molecule){
        return;
    }

    if(d->molecule){
        d->molecule->removeWatcher(this);
    }

    d->molecule = molecule;

    if(d->molecule){
        molecule->addWatcher(this);
    }
}

/// Returns the molecule that the watcher is monitoring.
const Molecule* MoleculeWatcher::molecule() const
{
    return d->molecule;
}

// --- Signals ------------------------------------------------------------- //
/// \fn void MoleculeWatcher::atomAdded(const Atom *atom)
///
/// This signal is emitted when an atom is added.

/// \fn void MoleculeWatcher::atomRemoved(const Atom *atom)
///
/// This signal is emitted when an atom is removed.

/// \fn void MoleculeWatcher::atomElementChanged(const Atom *atom)
///
/// This signal is emitted when an atom's element changes.

/// \fn void MoleculeWatcher::atomPositionChanged(const Atom *atom)
///
/// This signal is emitted when an atom's position changes.

/// \fn void MoleculeWatcher::bondAdded(const Bond *bond)
///
/// This signal is emitted when a bond is added.

/// \fn void MoleculeWatcher::bondRemoved(const Bond *bond)
///
/// This signal is emitted when a bond is removed.

/// \fn void MoleculeWatcher::bondOrderChanged(const Bond *bond)
///
/// This signal is emitted when a bond's order changes.

/// \fn void MoleculeWatcher::nameChanged(const Molecule *molecule)
///
/// This signal is emitted when the molecule's name changes.

// --- Events ---------------------------------------------------- //
void MoleculeWatcher::moleculeChanged(const Molecule *molecule, ChangeType changeType)
{
    switch(changeType){
        case NameChanged:
            nameChanged(molecule);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::atomChanged(const Atom *atom, ChangeType changeType)
{
    switch(changeType){
        case AtomAdded:
            atomAdded(atom);
            break;
        case AtomRemoved:
            atomRemoved(atom);
            break;
        case AtomElementChanged:
            atomElementChanged(atom);
            break;
        case AtomPositionChanged:
            atomPositionChanged(atom);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::bondChanged(const Bond *bond, ChangeType changeType)
{
    switch(changeType){
        case BondAdded:
            bondAdded(bond);
            break;
        case BondRemoved:
            bondRemoved(bond);
            break;
        case BondOrderChanged:
            bondOrderChanged(bond);
        default:
            break;
    }
}

} // end chemkit namespace
