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

#include "atom.h"
#include "foreach.h"
#include "vector3.h"
#include "geometry.h"
#include "molecule.h"
#include "cartesiancoordinates.h"
#include "coordinateset.h"

namespace chemkit {

// === MoleculeAlignerPrivate ============================================== //
class MoleculeAlignerPrivate
{
public:
    std::map<Atom *, Atom *> mapping;
    const Molecule *sourceMolecule;
    const Molecule *targetMolecule;
    const CoordinateSet *sourceCoordinates;
    const CoordinateSet *targetCoordinates;
};

// === MoleculeAligner ===================================================== //
/// \class MoleculeAligner moleculealigner.h chemkit/moleculealigner.h
/// \ingroup chemkit
/// \brief The MoleculeAligner class aligns two molecules based on
///        their atomic coordinates.
///
/// This class implements the \blueobeliskalgorithm{alignmentKabsch}.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new molecule aligner object using \p mapping.
MoleculeAligner::MoleculeAligner(const std::map<Atom *, Atom *> &mapping)
    : d(new MoleculeAlignerPrivate)
{
    setMapping(mapping);

    d->sourceCoordinates = 0;
    d->targetCoordinates = 0;
}

/// Create a new molecule aligner object using a mapping between the
/// \p source and \p target molecules.
MoleculeAligner::MoleculeAligner(const Molecule *source, const Molecule *target)
    : d(new MoleculeAlignerPrivate)
{
    d->sourceMolecule = source;
    d->targetMolecule = target;
    d->sourceCoordinates = 0;
    d->targetCoordinates = 0;

    int size = std::min(source->size(), target->size());

    for(int i = 0; i < size; i++){
        d->mapping[source->atom(i)] = target->atom(i);
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
    return d->sourceMolecule;
}

/// Returns the target molecule.
const Molecule* MoleculeAligner::targetMolecule() const
{
    return d->targetMolecule;
}

/// Sets the atom mapping to \p mapping.
void MoleculeAligner::setMapping(const std::map<Atom *, Atom *> &mapping)
{
    d->mapping = mapping;

    if(!mapping.empty()){
        Atom *a = mapping.begin()->first;
        Atom *b = mapping.begin()->second;

        d->sourceMolecule = a->molecule();
        d->targetMolecule = b->molecule();
    }
}

/// Returns the atom mapping.
std::map<Atom *, Atom *> MoleculeAligner::mapping() const
{
    return d->mapping;
}

/// Sets the coordinate set to use for the source molecule.
void MoleculeAligner::setSourceCoordinateSet(const CoordinateSet *coordinates)
{
    d->sourceCoordinates = coordinates;
}

/// Returns the coordinate set for the source molecule.
const CoordinateSet* MoleculeAligner::sourceCoordinateSet() const
{
    return d->sourceCoordinates;
}

/// Sets the coordinate set for the target molecule.
void MoleculeAligner::setTargetCoordinateSet(const CoordinateSet *coordinates)
{
    d->targetCoordinates = coordinates;
}

/// Returns the conformer used for the target molecule.
const CoordinateSet* MoleculeAligner::targetCoordinateSet() const
{
    return d->targetCoordinates;
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the root mean square deviation between the coordinates
/// of the source and target molecules.
Real MoleculeAligner::rmsd() const
{
    CartesianCoordinates *sourceMatrix = sourceCoordinates();
    CartesianCoordinates *targetMatrix = targetCoordinates();

    Real rmsd = this->rmsd(sourceMatrix, targetMatrix);

    delete sourceMatrix;
    delete targetMatrix;

    return rmsd;
}

/// Returns a 3x3 rotation matrix that represents the optimal
/// rotation of the source molecule to minimize the root mean square
/// deviation.
Eigen::Matrix<Real, 3, 3> MoleculeAligner::rotationMatrix() const
{
    CartesianCoordinates *sourceMatrix = sourceCoordinates();
    CartesianCoordinates *targetMatrix = targetCoordinates();

    sourceMatrix->moveBy(-sourceMatrix->center());
    targetMatrix->moveBy(-targetMatrix->center());

    Eigen::Matrix<Real, 3, 3> covarianceMatrix = targetMatrix->multiply(sourceMatrix);

    delete sourceMatrix;
    delete targetMatrix;

    Eigen::Matrix<Real, 3, 3> rotationMatrix = Eigen::Matrix<Real, 3, 3>::Identity();

    int d = covarianceMatrix.determinant() >= 0 ? 1 : -1;
    rotationMatrix(2, 2) = d;

    // compute singular value decomposition of the covariance matrix
    Eigen::JacobiSVD<Eigen::Matrix<Real, 3, 3> > svd(covarianceMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);

    rotationMatrix = svd.matrixU() * rotationMatrix * svd.matrixV().transpose();

    return rotationMatrix;
}

/// Returns a vector containing the displacement between the centers
/// of the source and target molecules.
Vector3 MoleculeAligner::displacementVector() const
{
    CartesianCoordinates *sourceMatrix = sourceCoordinates();
    CartesianCoordinates *targetMatrix = targetCoordinates();

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
    Eigen::Matrix<Real, 3, 3> matrix = rotationMatrix();
    foreach(Atom *atom, molecule->atoms()){
        atom->setPosition(matrix * atom->position());
    }

    Vector3 displacement = displacementVector();
    foreach(Atom *atom, molecule->atoms()){
        atom->setPosition(atom->position() + displacement);
    }
}

// --- Static Methods ------------------------------------------------------ //
/// Returns the root mean square deviation between the coordinates
/// in \p a and \p b.
Real MoleculeAligner::rmsd(const CartesianCoordinates *a, const CartesianCoordinates *b)
{
    int size = std::min(a->size(), b->size());

    Real sum = 0;

    for(int i = 0; i < size; i++){
        sum += chemkit::geometry::distanceSquared(a->position(i), b->position(i));
    }

    return sqrt(sum / size);
}

// --- Internal Methods ---------------------------------------------------- //
CartesianCoordinates* MoleculeAligner::sourceCoordinates() const
{
    if(d->sourceCoordinates){
        return new CartesianCoordinates(*d->sourceCoordinates->cartesianCoordinates());
    }
    else{
        return new CartesianCoordinates(*d->sourceMolecule->coordinates());
    }
}

CartesianCoordinates* MoleculeAligner::targetCoordinates() const
{
    if(d->targetCoordinates){
        return new CartesianCoordinates(*d->targetCoordinates->cartesianCoordinates());
    }
    else{
        return new CartesianCoordinates(*d->targetMolecule->coordinates());
    }
}

} // end chemkit namespace
