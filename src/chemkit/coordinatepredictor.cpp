/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "coordinatepredictor.h"

#include "atom.h"
#include "foreach.h"
#include "molecule.h"
#include "concurrent.h"

namespace chemkit {

// === CoordinatePredictorPrivate ========================================== //
class CoordinatePredictorPrivate
{
public:
    const Molecule *molecule;
};

// === CoordinatePredictor ================================================= //
/// \class CoordinatePredictor coordinatepredictor.h chemkit/coordinatepredictor.h
/// \ingroup chemkit
/// \brief The CoordinatePredictor class predicts 3D coordinates for
///        a molecule.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new coordinate predictor for \p molecule.
CoordinatePredictor::CoordinatePredictor(const Molecule *molecule)
    : d(new CoordinatePredictorPrivate)
{
    d->molecule = molecule;
}

/// Destroys the coordinate predictor object.
CoordinatePredictor::~CoordinatePredictor()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule to \p molecule.
void CoordinatePredictor::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the current molecule.
const Molecule* CoordinatePredictor::molecule() const
{
    return d->molecule;
}

// --- Static Methods ------------------------------------------------------ //
/// Predicts and assigns 3D coordinates for the atoms in \p molecule.
void CoordinatePredictor::predictCoordinates(Molecule *molecule)
{
    // store terminal hydrogens
    std::vector<Atom *> hydrogens;
    foreach(Atom *atom, molecule->atoms()){
        if(atom->isTerminalHydrogen()){
            hydrogens.push_back(atom);
        }
    }

    // assign heavy atom positions
    size_t heavyAtomCount = molecule->size() - hydrogens.size();

    foreach(Atom *atom, molecule->atoms()){
        if(atom->isTerminalHydrogen()){
            continue;
        }

        atom->setPosition(heavyAtomCount *
                          chemkit::Vector3::Random().normalized());
    }

    // assign hydrogen positions
    foreach(Atom *hydrogen, hydrogens){
        const Atom *neighbor = hydrogen->neighbor(0);

        hydrogen->setPosition(neighbor->position() +
                              chemkit::Vector3::Random().normalized());
    }
}

/// Runs the predictCoordinates() method asynchronously and returns
/// a future containing the result.
///
/// \internal
boost::shared_future<void> CoordinatePredictor::predictCoordinatesAsync(Molecule *molecule)
{
  return chemkit::concurrent::run(
      boost::bind(&CoordinatePredictor::predictCoordinates, molecule));
}

/// Adjusts the coordinates of the atoms in \p molecule to ensure that
/// no two atoms are within \p distance Angstroms of each other. Returns
/// \c true if at least one close contact was found and eliminated.
bool CoordinatePredictor::eliminateCloseContacts(Molecule *molecule, Real distance)
{
    bool done = false;
    bool modified = false;

    while(!done){
        done = true;

        for(size_t i = 0; i < molecule->size(); i++){
            Atom *a = molecule->atom(i);

            for(size_t j = i + 1; j < molecule->size(); j++){
                Atom *b = molecule->atom(j);

                if(a->distance(b) < distance){
                    done = false;

                    // move atom b by a random unit vector
                    b->setPosition(b->position() +
                                   distance * Vector3::Random().normalized());

                    // set modified flag
                    modified = true;
                }
            }
        }
    }

    return modified;
}

} // end chemkit namespace
