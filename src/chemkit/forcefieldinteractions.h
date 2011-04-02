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

#ifndef CHEMKIT_FORCEFIELDINTERACTIONS_H
#define CHEMKIT_FORCEFIELDINTERACTIONS_H

#include "chemkit.h"

#include <vector>

namespace chemkit {

class Atom;
class Molecule;
class ForceField;
class ForceFieldAtom;
class ForceFieldInteractionsPrivate;

class CHEMKIT_EXPORT ForceFieldInteractions
{
    public:
        // construction and destruction
        ForceFieldInteractions(const Molecule *molecule, const ForceField *forceField);
        ~ForceFieldInteractions();

        // properties
        const Molecule* molecule() const;
        const ForceField* forceField() const;

        // interactions
        std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > bondedPairs();
        std::vector<std::vector<const ForceFieldAtom *> > angleGroups();
        std::vector<std::vector<const ForceFieldAtom *> > torsionGroups();
        std::vector<std::pair<const ForceFieldAtom *, const ForceFieldAtom *> > nonbondedPairs();

    private:
        bool atomsWithinTwoBonds(const Atom *a, const Atom *b);

    private:
        ForceFieldInteractionsPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDINTERACTIONS_H
