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

#include "moleculealigner.h"

#include "vector3.h"
#include "molecule.h"
#include "coordinates.h"

namespace chemkit {

// === MoleculeAlignerPrivate ============================================== //
class MoleculeAlignerPrivate
{
    public:
        int size;
        AtomMapping mapping;
        const Conformer *sourceConformer;
        const Conformer *targetConformer;
};

// === MoleculeAligner ===================================================== //
/// \class MoleculeAligner moleculealigner.h chemkit/moleculealigner.h
/// \ingroup chemkit
/// \brief The MoleculeAligner class aligns two molecules based on
///        their atomic coordinates.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new molecule aligner object using \p mapping.
MoleculeAligner::MoleculeAligner(const AtomMapping &mapping)
    : d(new MoleculeAlignerPrivate)
{
    d->mapping = mapping;
    d->sourceConformer = 0;
    d->targetConformer = 0;
    d->size = mapping.size();
}

/// Create a new molecule aligner object using a mapping between the
/// \p source and \p target molecules.
MoleculeAligner::MoleculeAligner(const Molecule *source, const Molecule *target)
    : d(new MoleculeAlignerPrivate)
{
    d->mapping = AtomMapping(source, target);
    d->sourceConformer = 0;
    d->targetConformer = 0;
    d->size = qMin(source->size(), target->size());

    for(int i = 0; i < d->size; i++){
        d->mapping.add(source->atom(i), target->atom(i));
    }
}

/// Destroys the molecule aligner object.
MoleculeAligner::~MoleculeAligner()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the source molecule.
const Molecule* MoleculeAligner::sourceMolecule() const
{
    return d->mapping.source();
}

/// Returns the target molecule.
const Molecule* MoleculeAligner::targetMolecule() const
{
    return d->mapping.target();
}

/// Sets the atom mapping to \p mapping.
void MoleculeAligner::setMapping(const AtomMapping &mapping)
{
    d->mapping = mapping;
}

/// Returns the atom mapping.
const AtomMapping& MoleculeAligner::mapping() const
{
    return d->mapping;
}

/// Sets the conformer to use for the source molecule.
void MoleculeAligner::setSourceConformer(const Conformer *conformer)
{
    if(conformer->molecule() != sourceMolecule())
        return;

    d->sourceConformer = conformer;
}

/// Returns the conformer for the source molecule.
const Conformer* MoleculeAligner::sourceConformer() const
{
    return d->sourceConformer;
}

/// Sets the conformer for the target molecule.
void MoleculeAligner::setTargetConformer(const Conformer *conformer)
{
    if(conformer->molecule() != targetMolecule())
        return;

    d->targetConformer = conformer;
}

/// Returns the conformer used for the target molecule.
const Conformer* MoleculeAligner::targetConformer() const
{
    return d->targetConformer;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the root mean square deviation between the coordinates
/// of the source and target molecules.
Float MoleculeAligner::deviation() const
{
    Coordinates *sourceMatrix = sourceCoordinates();
    Coordinates *targetMatrix = targetCoordinates();

    Float rmsd = this->rmsd(sourceMatrix, targetMatrix);

    delete sourceMatrix;
    delete targetMatrix;

    return rmsd;
}

/// Returns a 3x3 rotation matrix that represents the optimal
/// rotation of the source molecule to minimize the root mean square
/// deviation.
StaticMatrix<Float, 3, 3> MoleculeAligner::rotationMatrix() const
{
    Coordinates *sourceMatrix = sourceCoordinates();
    Coordinates *targetMatrix = targetCoordinates();

    sourceMatrix->moveBy(-sourceMatrix->center());
    targetMatrix->moveBy(-targetMatrix->center());

    StaticMatrix<Float, 3, 3> covarianceMatrix = targetMatrix->multiply(sourceMatrix);

    delete sourceMatrix;
    delete targetMatrix;

    StaticMatrix<Float, 3, 3> rotationMatrix = StaticMatrix<Float, 3, 3>::identity();

    int d = covarianceMatrix.determinant() >= 0 ? 1 : -1;
    rotationMatrix(2, 2) = d;

    // compute singular value decomposition of the covariance matrix
    StaticMatrix<Float, 3, 3> U;
    StaticMatrix<Float, 3, 3> Vt;
    StaticVector<Float, 3> S;

    covarianceMatrix.svd(&U, &S, &Vt);

    rotationMatrix = U * rotationMatrix * Vt;

    return rotationMatrix;
}

/// Returns a vector containing the displacement between the centers
/// of the source and target molecules.
Vector3 MoleculeAligner::displacementVector() const
{
    Coordinates *sourceMatrix = sourceCoordinates();
    Coordinates *targetMatrix = targetCoordinates();

    Vector3 displacement = targetCoordinates()->center() - sourceCoordinates()->center();

    delete sourceMatrix;
    delete targetMatrix;

    return displacement;
}

/// Aligns the molecule by transforming it by the rotation matrix
/// returned rotationMatrix() and moving it by the vector returned
/// from displacementVector().
void MoleculeAligner::align(Molecule *molecule)
{
    StaticMatrix<Float, 3, 3> matrix = rotationMatrix();
    foreach(Atom *atom, molecule->atoms()){
        atom->setPosition(matrix.multiply(atom->position()));
    }

    Vector3 displacement = displacementVector();
    foreach(Atom *atom, molecule->atoms()){
        atom->moveBy(displacement);
    }
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the root mean square deviation between the coordinates
/// in \p a and \p b.
Float MoleculeAligner::rmsd(const Coordinates *a, const Coordinates *b)
{
    int size = qMin(a->size(), b->size());

    Float sum = 0;

    for(int i = 0; i < size; i++){
        sum += Point3::distanceSquared(a->position(i), b->position(i));
    }

    return sqrt(sum / size);
}

// --- Internal Methods ---------------------------------------------------- //
Coordinates* MoleculeAligner::sourceCoordinates() const
{
    if(d->sourceConformer){
        return new Coordinates(d->sourceConformer);
    }
    else{
        return new Coordinates(sourceMolecule());
    }
}

Coordinates* MoleculeAligner::targetCoordinates() const
{
    if(d->targetConformer){
        return new Coordinates(d->targetConformer);
    }
    else{
        return new Coordinates(targetMolecule());
    }
}

} // end chemkit namespace
