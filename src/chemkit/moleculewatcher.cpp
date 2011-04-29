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
    : QObject(),
      d(new MoleculeWatcherPrivate)
{
    d->molecule = molecule;

    if(molecule){
        molecule->addWatcher(this);
    }
}

/// Destroys the molecule watcher object.
MoleculeWatcher::~MoleculeWatcher()
{
    if(d->molecule){
        d->molecule->removeWatcher(this);
    }

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule to be watched.
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
        d->molecule->addWatcher(this);
    }
}

/// Returns the molecule that is being watched.
const Molecule* MoleculeWatcher::molecule() const
{
    return d->molecule;
}

// --- Signals ------------------------------------------------------------- //
/// \fn void MoleculeWatcher::atomAdded(const chemkit::Atom *atom)
///
/// This signal is emitted when an atom is added.

/// \fn void MoleculeWatcher::atomRemoved(const chemkit::Atom *atom)
///
/// This signal is emitted when an atom is removed.

/// \fn void MoleculeWatcher::atomAtomicNumberChanged(const chemkit::Atom *atom)
///
/// This signal is emitted when an atom's atomic number changes.

/// \fn void MoleculeWatcher::atomPositionChanged(const chemkit::Atom *atom)
///
/// This signal is emitted when an atom's position changes.

/// \fn void MoleculeWatcher::bondAdded(const chemkit::Bond *bond)
///
/// This signal is emitted when a bond is added.

/// \fn void MoleculeWatcher::bondRemoved(const chemkit::Bond *bond)
///
/// This signal is emitted when a bond is removed.

/// \fn void MoleculeWatcher::bondOrderChanged(const chemkit::Bond *bond)
///
/// This signal is emitted when a bond's order changes.

/// \fn void MoleculeWatcher::residueAdded(const chemkit::Residue *residue)
///
/// This signal is emitted when a residue is added.

/// \fn void MoleculeWatcher::residueRemoved(const chemkit::Residue *residue)
///
/// This signal is emitted when a residue is removed.

/// \fn void MoleculeWatcher::conformerAdded(const chemkit::Conformer *conformer)
///
/// This signal is emitted when a conformer is added.

/// \fn void MoleculeWatcher::conformerRemoved(const chemkit::Conformer *conformer)
///
/// This signal is emitted when a conformer is removed.

/// \fn void MoleculeWatcher::nameChanged(const chemkit::Molecule *molecule)
///
/// This signal is emitted when the molecule's name changes.

// --- Internal Methods ---------------------------------------------------- //
void MoleculeWatcher::notifyObservers(const Molecule *molecule, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::NameChanged:
            Q_EMIT nameChanged(molecule);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Atom *atom, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::AtomAdded:
            Q_EMIT atomAdded(atom);
            break;
        case Molecule::AtomRemoved:
            Q_EMIT atomRemoved(atom);
            break;
        case Molecule::AtomAtomicNumberChanged:
            Q_EMIT atomAtomicNumberChanged(atom);
            break;
        case Molecule::AtomPositionChanged:
            Q_EMIT atomPositionChanged(atom);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Bond *bond, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::BondAdded:
            Q_EMIT bondAdded(bond);
            break;
        case Molecule::BondRemoved:
            Q_EMIT bondRemoved(bond);
            break;
        case Molecule::BondOrderChanged:
            Q_EMIT bondOrderChanged(bond);
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Conformer *conformer, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::ConformerAdded:
            Q_EMIT conformerAdded(conformer);
            break;
        case Molecule::ConformerRemoved:
            Q_EMIT conformerRemoved(conformer);
            break;
        default:
            break;
    }
}

} // end chemkit namespace
