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

#include "partialchargemodel.h"

#include "atom.h"
#include "foreach.h"
#include "molecule.h"
#include "pluginmanager.h"

namespace chemkit {

// === PartialChargeModelPrivate =========================================== //
class PartialChargeModelPrivate
{
public:
    std::string name;
    const Molecule *molecule;
};

// === PartialChargeModel ================================================== //
/// \class PartialChargeModel partialchargemodel.h chemkit/partialchargemodel.h
/// \ingroup chemkit
/// \brief The PartialChargeModel class provides a generic interface
///        to various partial charge models.
///
/// A list of supported partial charge models is available at:
/// http://wiki.chemkit.org/Features#Partial_Charge_Models

// --- Construction and Destruction ---------------------------------------- //
PartialChargeModel::PartialChargeModel(const std::string &name)
    : d(new PartialChargeModelPrivate)
{
    d->name = name;
    d->molecule = 0;
}

/// Destroys the partial charge model object.
PartialChargeModel::~PartialChargeModel()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the partial charge model.
std::string PartialChargeModel::name() const
{
    return d->name;
}

/// Sets the molecule for the model to \p molecule.
void PartialChargeModel::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the molecule for the model.
const Molecule* PartialChargeModel::molecule() const
{
    return d->molecule;
}

// --- Partial Charges ----------------------------------------------------- //
/// Returns the partial charge for \p atom.
Real PartialChargeModel::partialCharge(const Atom *atom) const
{
    CHEMKIT_UNUSED(atom);

    return 0;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new partial charge model with \p name. Returns \c 0
/// if \p name is invalid.
PartialChargeModel* PartialChargeModel::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<PartialChargeModel>(name);
}

/// Returns a list of available partial charge models.
std::vector<std::string> PartialChargeModel::models()
{
    return PluginManager::instance()->pluginClassNames<PartialChargeModel>();
}

/// This static convenience function assigns partial charges for atoms
/// in molecule using the specified \p model.
///
/// For example, to assign gasteiger partial charges for a molecule:
/// \code
/// PartialChargeModel::assignPartialCharges(molecule, "gasteiger");
/// \endcode
bool PartialChargeModel::assignPartialCharges(Molecule *molecule, const std::string &model)
{
    boost::scoped_ptr<PartialChargeModel> partialChargeModel(create(model));
    if(!partialChargeModel){
        return false;
    }

    partialChargeModel->setMolecule(molecule);

    foreach(Atom *atom, molecule->atoms()){
        atom->setPartialCharge(partialChargeModel->partialCharge(atom));
    }

    return true;
}

} // end chemkit namespace
