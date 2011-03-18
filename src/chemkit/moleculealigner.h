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

#ifndef CHEMKIT_MOLECULEALIGNER_H
#define CHEMKIT_MOLECULEALIGNER_H

#include "chemkit.h"

#include "vector3.h"
#include "atommapping.h"
#include "staticmatrix.h"

namespace chemkit {

class Molecule;
class Conformer;
class Coordinates;
class MoleculeAlignerPrivate;

class CHEMKIT_EXPORT MoleculeAligner
{
    public:
        // construction and destruction
        MoleculeAligner(const AtomMapping &mapping);
        MoleculeAligner(const Molecule *source, const Molecule *target);
        ~MoleculeAligner();

        // properties
        const Molecule* sourceMolecule() const;
        const Molecule* targetMolecule() const;
        void setMapping(const AtomMapping &mapping);
        const AtomMapping& mapping() const;
        void setSourceConformer(const Conformer *conformer);
        const Conformer* sourceConformer() const;
        void setTargetConformer(const Conformer *conformer);
        const Conformer* targetConformer() const;

        // geometry
        Float deviation() const;
        StaticMatrix<Float, 3, 3> rotationMatrix() const;
        Vector3 displacementVector() const;
        void align(Molecule *molecule);

        // static methods
        static Float rmsd(const Coordinates *a, const Coordinates *b);

    private:
        Coordinates *sourceCoordinates() const;
        Coordinates *targetCoordinates() const;

    private:
        MoleculeAlignerPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEALIGNER_H
