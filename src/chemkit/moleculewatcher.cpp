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
            emit nameChanged(molecule);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Atom *atom, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::AtomAdded:
            emit atomAdded(atom);
            break;
        case Molecule::AtomRemoved:
            emit atomRemoved(atom);
            break;
        case Molecule::AtomAtomicNumberChanged:
            emit atomAtomicNumberChanged(atom);
            break;
        case Molecule::AtomPositionChanged:
            emit atomPositionChanged(atom);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Bond *bond, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::BondAdded:
            emit bondAdded(bond);
            break;
        case Molecule::BondRemoved:
            emit bondRemoved(bond);
            break;
        case Molecule::BondOrderChanged:
            emit bondOrderChanged(bond);
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Residue *residue, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::ResidueAdded:
            emit residueAdded(residue);
            break;
        case Molecule::ResidueRemoved:
            emit residueRemoved(residue);
            break;
        default:
            break;
    }
}

void MoleculeWatcher::notifyObservers(const Conformer *conformer, Molecule::ChangeType changeType)
{
    switch(changeType){
        case Molecule::ConformerAdded:
            emit conformerAdded(conformer);
            break;
        case Molecule::ConformerRemoved:
            emit conformerRemoved(conformer);
            break;
        default:
            break;
    }
}

} // end chemkit namespace
