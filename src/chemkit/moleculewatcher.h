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

#ifndef CHEMKIT_MOLECULEWATCHER_H
#define CHEMKIT_MOLECULEWATCHER_H

#include "chemkit.h"

#include "molecule.h"

namespace chemkit {

class Atom;
class Bond;
class Residue;
class Molecule;
class Conformer;
class MoleculeWatcherPrivate;

class CHEMKIT_EXPORT MoleculeWatcher : public QObject
{
    Q_OBJECT

    public:
        // construction and destruction
        MoleculeWatcher(const Molecule *molecule = 0);
        ~MoleculeWatcher();

        // properties
        void setMolecule(const Molecule *molecule);
        const Molecule* molecule() const;

    signals:
        void atomAdded(const chemkit::Atom *atom);
        void atomRemoved(const chemkit::Atom *atom);
        void atomAtomicNumberChanged(const chemkit::Atom *atom);
        void atomPositionChanged(const chemkit::Atom *atom);
        void bondAdded(const chemkit::Bond *bond);
        void bondRemoved(const chemkit::Bond *bond);
        void bondOrderChanged(const chemkit::Bond *bond);
        void residueAdded(const chemkit::Residue *residue);
        void residueRemoved(const chemkit::Residue *residue);
        void conformerAdded(const chemkit::Conformer *conformer);
        void conformerRemoved(const chemkit::Conformer *conformer);
        void nameChanged(const chemkit::Molecule *molecule);

    private:
        void notifyObservers(const Molecule *molecule, Molecule::ChangeType changeType);
        void notifyObservers(const Atom *atom, Molecule::ChangeType changeType);
        void notifyObservers(const Bond *bond, Molecule::ChangeType changeType);
        void notifyObservers(const Residue *residue, Molecule::ChangeType changeType);
        void notifyObservers(const Conformer *conformer, Molecule::ChangeType changeType);

        friend class Molecule;

    private:
        MoleculeWatcherPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEWATCHER_H
