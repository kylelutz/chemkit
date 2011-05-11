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

#include "moleculealigner.h"

#include "vector3.h"
#include "geometry.h"
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
    d->size = std::min(source->size(), target->size());

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
    Q_FOREACH(Atom *atom, molecule->atoms()){
        atom->setPosition(matrix.multiply(atom->position()));
    }

    Vector3 displacement = displacementVector();
    Q_FOREACH(Atom *atom, molecule->atoms()){
        atom->moveBy(displacement);
    }
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the root mean square deviation between the coordinates
/// in \p a and \p b.
Float MoleculeAligner::rmsd(const Coordinates *a, const Coordinates *b)
{
    int size = std::min(a->size(), b->size());

    Float sum = 0;

    for(int i = 0; i < size; i++){
        sum += chemkit::geometry::distanceSquared(a->position(i), b->position(i));
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
