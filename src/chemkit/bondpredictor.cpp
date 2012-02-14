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

#include "bondpredictor.h"

#include "atom.h"
#include "foreach.h"
#include "molecule.h"

namespace chemkit {

// === BondPredictorPrivate ================================================ //
class BondPredictorPrivate
{
public:
    Molecule *molecule;
    Real tolerance;
    Real minimumBondLength;
    Real maximumBondLength;
};

// === BondPredictor ======================================================= //
/// \class BondPredictor bondpredictor.h chemkit/bondpredictor.h
/// \ingroup chemkit
/// \brief The BondPredictor class predicts bonds in a molecule.
///
/// The BondPredictor class predicts bonds for a molecule based on
/// the 3D coordinates of its atoms.
///
/// The easiest way to predict bonds for a molecule is by using the
/// static predictBonds() method as the following example shows:
/// \code
/// Molecule *molecule = ...
///
/// BondPredictor::predictBonds(molecule);
/// \endcode
///
/// This class implements the \blueobeliskalgorithm{rebondFrom3DCoordinates}.

/// \typedef BondPredictor::PredictedBond;
/// This tuple contains information about each predicted bond.
///
/// For example, the following code will retrieve each atom and
/// the bond order for the predicted bond:
/// \code
/// BondPredictor::PredictedBond bond = bondPredictor.predictedBonds()[0];
///
/// Atom *a = boost::get<0>(bond);
/// Atom *b = boost::get<1>(bond);
/// Bond::BondOrderType order = boost::get<2>(bond);
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
/// Create a new bond predictor object for molecule.
BondPredictor::BondPredictor(Molecule *molecule)
    : d(new BondPredictorPrivate)
{
    d->molecule = molecule;

    // default parameters
    d->minimumBondLength = 0.4;
    d->maximumBondLength = 5.0;
    d->tolerance = 0.45;
}

/// Destroys the bond predictor object.
BondPredictor::~BondPredictor()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the bond distance tolerance.
void BondPredictor::setTolerance(Real tolerance)
{
    d->tolerance = tolerance;
}

/// Returns the bond distance tolerance.
Real BondPredictor::tolerance() const
{
    return d->tolerance;
}

/// Sets the minimum bond length.
void BondPredictor::setMinimumBondLength(Real length)
{
    d->minimumBondLength = length;
}

/// Returns the minimum bond length.
Real BondPredictor::minimumBondLength() const
{
    return d->minimumBondLength;
}

/// Sets the maximum bond length.
void BondPredictor::setMaximumBondLength(Real length)
{
    d->maximumBondLength = length;
}

/// Returns the maximum bond length.
Real BondPredictor::maximumBondLength() const
{
    return d->maximumBondLength;
}

/// Returns the molecule.
Molecule* BondPredictor::molecule() const
{
    return d->molecule;
}

// --- Prediction ---------------------------------------------------------- //
/// Returns a list of pairs of atoms that are predicted to be bonded.
std::vector<BondPredictor::PredictedBond> BondPredictor::predictedBonds()
{
    std::vector<PredictedBond> bonds;

    if(!d->molecule)
        return bonds;

    std::vector<Atom *> atoms(d->molecule->atoms().begin(), d->molecule->atoms().end());

    for(unsigned int i = 0; i < atoms.size(); i++){
        for(unsigned int j = i+1; j < atoms.size(); j++){
            if(couldBeBonded(atoms[i], atoms[j])){
                bonds.push_back(boost::make_tuple(atoms[i], atoms[j], Bond::Single));
            }
        }
    }

    return bonds;
}

// --- Static Methods ------------------------------------------------------ //
/// This static convenience method predicts the bonds for \p molecule
/// and adds each predicted bond with the Molecule::addBond() method.
void BondPredictor::predictBonds(Molecule *molecule)
{
    BondPredictor predictor(molecule);

    foreach(const PredictedBond &bond, predictor.predictedBonds()){
        molecule->addBond(boost::get<0>(bond), boost::get<1>(bond), boost::get<2>(bond));
    }
}

// --- Internal Methods ---------------------------------------------------- //
// Returns \c true if the atoms could feasibly be bonded.
bool BondPredictor::couldBeBonded(Atom *a, Atom *b) const
{
    Real distance = a->distance(b);

    if(distance > minimumBondLength() &&
       distance < maximumBondLength() &&
       std::abs((a->covalentRadius() + b->covalentRadius()) - distance) < tolerance())
        return true;
    else
        return false;
}

} // end chemkit namespace
