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

#include "partialchargepredictor.h"

#include "atom.h"
#include "molecule.h"
#include "pluginmanager.h"

namespace chemkit {

namespace {

QHash<QString, PartialChargePredictor::CreateFunction> predictorPlugins;

} // end anonymous namespace

// === PartialChargePredictorPrivate ======================================= //
class PartialChargePredictorPrivate
{
    public:
        QString name;
        const Molecule *molecule;
};

// === PartialChargePredictor ============================================== //
/// \class PartialChargePredictor partialchargepredictor.h chemkit/partialchargepredictor.h
/// \ingroup chemkit
/// \brief The PartialChargePredictor class provides a generic
///        interface to partial charge prediction algorithms.

// --- Construction and Destruction ---------------------------------------- //
PartialChargePredictor::PartialChargePredictor(const QString &name)
    : d(new PartialChargePredictorPrivate)
{
    d->name = name;
    d->molecule = 0;
}

/// Destroys the partial charge predictor object.
PartialChargePredictor::~PartialChargePredictor()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the partial charge predictor.
QString PartialChargePredictor::name() const
{
    return d->name;
}

/// Sets the molecule for the predictor to \p molecule.
void PartialChargePredictor::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;

    assignPartialCharges(molecule);
}

/// Returns the molecule for the predictor.
const Molecule* PartialChargePredictor::molecule() const
{
    return d->molecule;
}

// --- Partial Charges ----------------------------------------------------- //
/// Returns the partial charge for the atom at \p index.
Float PartialChargePredictor::partialCharge(int index) const
{
    Q_UNUSED(index);

    return 0;
}

/// Returns the partial charge for \p atom.
Float PartialChargePredictor::partialCharge(const Atom *atom) const
{
    return partialCharge(atom->index());
}

void PartialChargePredictor::assignPartialCharges(const Molecule *molecule)
{
    Q_UNUSED(molecule);
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new partial charge predictor with \p name. Returns \c 0
/// if \p name is invalid.
PartialChargePredictor* PartialChargePredictor::create(const QString &name)
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    CreateFunction createFunction = predictorPlugins.value(name.toLower());
    if(createFunction){
        return createFunction();
    }

    return 0;
}

/// Returns a list of available partial charge predictors.
QStringList PartialChargePredictor::predictors()
{
    // ensure default plugins are loaded
    PluginManager::instance()->loadDefaultPlugins();

    return predictorPlugins.keys();
}

void PartialChargePredictor::registerPredictor(const QString &name, CreateFunction function)
{
    predictorPlugins.insert(name.toLower(), function);
}

void PartialChargePredictor::unregisterPredictor(const QString &name, CreateFunction function)
{
    if(predictorPlugins.value(name.toLower()) == function){
        predictorPlugins.remove(name.toLower());
    }
}

bool PartialChargePredictor::predictPartialCharges(Molecule *molecule, const QString &predictorName)
{
    PartialChargePredictor *predictor = create(predictorName);
    if(!predictor){
        return false;
    }

    predictor->setMolecule(molecule);

    foreach(chemkit::Atom *atom, molecule->atoms()){
        atom->setPartialCharge(predictor->partialCharge(atom));
    }

    delete predictor;

    return true;
}

} // end chemkit namespace
