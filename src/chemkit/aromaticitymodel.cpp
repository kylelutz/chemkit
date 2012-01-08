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

#include "aromaticitymodel.h"

#include "atom.h"
#include "bond.h"
#include "foreach.h"
#include "molecule.h"
#include "pluginmanager.h"

namespace chemkit {

// === AromaticityModelPrivate ============================================= //
class AromaticityModelPrivate
{
public:
    std::string name;
    const Molecule *molecule;
};

// === AromaticityModel ==================================================== //
/// \class AromaticityModel aromaticitymodel.h chemkit/aromaticitymodel.h
/// \ingroup chemkit
/// \brief The AromaticityModel class represents a model of
///        aromaticity.
///
/// A list of supported aromaticity models is available at:
/// http://wiki.chemkit.org/Features#Aromaticity_Models

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new aromaticity model.
AromaticityModel::AromaticityModel()
    : d(new AromaticityModelPrivate)
{
    d->molecule = 0;
}

/// Creates a new aromaticity model with \p name.
AromaticityModel::AromaticityModel(const std::string &name)
    : d(new AromaticityModelPrivate)
{
    d->name = name;
    d->molecule = 0;
}

/// Destroys the aromaticity model object.
AromaticityModel::~AromaticityModel()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the aromaticity model.
std::string AromaticityModel::name() const
{
    return d->name;
}

/// Sets the molecule to \p molecule.
void AromaticityModel::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the current molecule.
const Molecule* AromaticityModel::molecule() const
{
    return d->molecule;
}

// --- Aromaticity --------------------------------------------------------- //
/// Returns \c true if \p atom is aromatic according to the model.
bool AromaticityModel::isAromatic(const Atom *atom) const
{
    return isAromaticAtom(atom);
}

/// Returns \c true if \p bond is aromatic according to the model.
bool AromaticityModel::isAromatic(const Bond *bond) const
{
    return isAromaticBond(bond);
}

/// Returns \c true if \p ring is aromatic according to the model.
bool AromaticityModel::isAromatic(const Ring *ring) const
{
    return isAromaticRing(ring);
}

/// Returns \c true if \p atom is aromatic. This method can be
/// reimplemented by each aromaticity model.
///
/// The default implementation returns \c true if the atom is a
/// member of any aromatic ring as determined by isAromaticRing().
bool AromaticityModel::isAromaticAtom(const Atom *atom) const
{
    foreach(const Ring *ring, atom->rings()){
        if(isAromaticRing(ring)){
            return true;
        }
    }

    return false;
}

/// Returns \c true if \p bond is aromatic. This method can be
/// reimplemented by each aromaticity model.
///
/// The default implementation returns \c true if the atom is a
/// member of any aromatic ring as determined by isAromaticRing().
bool AromaticityModel::isAromaticBond(const Bond *bond) const
{
    foreach(const Ring *ring, bond->rings()){
        if(isAromaticRing(ring)){
            return true;
        }
    }

    return false;
}

/// Returns \c true if \p ring is aromatic. This method should be
/// reimplemented by each aromaticity model.
///
/// The default implementation returns \c false.
bool AromaticityModel::isAromaticRing(const Ring *ring) const
{
    CHEMKIT_UNUSED(ring);

    return false;
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new aromaticity model object. Returns \c 0 if \p name
/// is invalid.
AromaticityModel* AromaticityModel::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<AromaticityModel>(name);
}

/// Returns a list of the supported aromaticity models.
std::vector<std::string> AromaticityModel::models()
{
    return PluginManager::instance()->pluginClassNames<AromaticityModel>();
}

} // end chemkit namespace
