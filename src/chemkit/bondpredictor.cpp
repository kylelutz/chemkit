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

#include "bondpredictor.h"

#include "atom.h"
#include "molecule.h"

namespace chemkit {

// === BondPredictorPrivate ================================================ //
class BondPredictorPrivate
{
    public:
        Molecule *molecule;
        Float tolerance;
        Float minimumBondLength;
        Float maximumBondLength;
};

// === BondPredictor ======================================================= //
/// \class BondPredictor bondpredictor.h chemkit/bondpredictor.h
/// \ingroup chemkit
/// \brief The BondPredictor class predicts bonds in a molecule.

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
void BondPredictor::setTolerance(Float tolerance)
{
    d->tolerance = tolerance;
}

/// Returns the bond distance tolerance.
Float BondPredictor::tolerance() const
{
    return d->tolerance;
}

/// Sets the minimum bond length.
void BondPredictor::setMinimumBondLength(Float length)
{
    d->minimumBondLength = length;
}

/// Returns the minimum bond length.
Float BondPredictor::minimumBondLength() const
{
    return d->minimumBondLength;
}

/// Sets the maximum bond length.
void BondPredictor::setMaximumBondLength(Float length)
{
    d->maximumBondLength = length;
}

/// Returns the maximum bond length.
Float BondPredictor::maximumBondLength() const
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
QList<std::pair<Atom *, Atom *> > BondPredictor::predictedBonds()
{
    QList<std::pair<Atom *, Atom *> >  bonds;

    if(!d->molecule)
        return bonds;

    QList<Atom *> atoms = d->molecule->atoms();

    for(int i = 0; i < atoms.size(); i++){
        for(int j = i+1; j < atoms.size(); j++){
            if(couldBeBonded(atoms[i], atoms[j])){
                bonds.append(std::make_pair(atoms[i], atoms[j]));
            }
        }
    }

    return bonds;
}

/// Returns \c true if the atoms could feasibly be bonded.
bool BondPredictor::couldBeBonded(Atom *a, Atom *b) const
{
    Float distance = a->distance(b);

    if(distance > minimumBondLength() &&
       distance < maximumBondLength() &&
       qAbs((a->covalentRadius() + b->covalentRadius()) - distance) < tolerance())
        return true;
    else
        return false;
}

// --- Static Methods ------------------------------------------------------ //
/// Predict bonds for the molecule.
void BondPredictor::predictBonds(Molecule *molecule)
{
    BondPredictor predictor(molecule);

    std::pair<Atom *, Atom *> bond;
    Q_FOREACH(bond, predictor.predictedBonds()){
        molecule->addBond(bond.first, bond.second);
    }
}

} // end chemkit namespace
