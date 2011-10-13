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

#ifndef CHEMKIT_MOLECULEALIGNER_H
#define CHEMKIT_MOLECULEALIGNER_H

#include "chemkit.h"

#include <map>

#include <Eigen/Core>

#include "vector3.h"

namespace chemkit {

class Atom;
class Molecule;
class Conformer;
class Coordinates;
class MoleculeAlignerPrivate;

class CHEMKIT_EXPORT MoleculeAligner
{
    public:
        // construction and destruction
        MoleculeAligner(const std::map<Atom *, Atom *> &mapping);
        MoleculeAligner(const Molecule *source, const Molecule *target);
        ~MoleculeAligner();

        // properties
        const Molecule* sourceMolecule() const;
        const Molecule* targetMolecule() const;
        void setMapping(const std::map<Atom *, Atom *> &mapping);
        std::map<Atom *, Atom *> mapping() const;
        void setSourceConformer(const Conformer *conformer);
        const Conformer* sourceConformer() const;
        void setTargetConformer(const Conformer *conformer);
        const Conformer* targetConformer() const;

        // geometry
        Real deviation() const;
        Eigen::Matrix<Real, 3, 3> rotationMatrix() const;
        Vector3 displacementVector() const;
        void align(Molecule *molecule);

        // static methods
        static Real rmsd(const Coordinates *a, const Coordinates *b);

    private:
        Coordinates *sourceCoordinates() const;
        Coordinates *targetCoordinates() const;

    private:
        MoleculeAlignerPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEALIGNER_H
