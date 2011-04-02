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

#include "moleculardescriptor.h"

#include <map>
#include <boost/algorithm/string.hpp>

#include "foreach.h"
#include "pluginmanager.h"

namespace chemkit {

// === MolecularDescriptorPrivate ========================================== //
class MolecularDescriptorPrivate
{
    public:
        std::string name;
};

// === MolecularDescriptor ================================================= //
/// \class MolecularDescriptor moleculardescriptor.h chemkit/moleculardescriptor.h
/// \ingroup chemkit
/// \brief The MolecularDescriptor class provides a generic interface
///        for the calculation of molecular descriptors.
///
/// The following molecular descriptors are supported in chemkit:
///     - \c randic-index
///     - \c wiener-index
///     - \c rotatable-bonds
///     - \c hydrogen-bond-donors
///     - \c hydrogen-bond-acceptors
///
/// For example, to calcuate the wiener index of a molecule:
/// \code
/// // load molecule from file, string, etc.
/// Molecule *molecule = ...;
///
/// // create wiener index descriptor
/// MolecularDescriptor *descriptor = MolecularDescriptor::create("wiener-index");
/// if(!descriptor){
///     // wiener index descriptor not available
///     return;
/// }
///
/// // calculate wiener index
/// int wienerIndex = descriptor->value(molecule).toInt();
///
/// // destroy descriptor object
/// delete descriptor;
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecular descriptor object.
MolecularDescriptor::MolecularDescriptor(const std::string &name)
    : d(new MolecularDescriptorPrivate)
{
    d->name = name;
}

/// Destroys the molecular descriptor object.
MolecularDescriptor::~MolecularDescriptor()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the descriptor.
std::string MolecularDescriptor::name() const
{
    return d->name;
}

// --- Descriptor ---------------------------------------------------------- //
/// Calculates the value of the descriptor for \p molecule.
QVariant MolecularDescriptor::value(const Molecule *molecule) const
{
    Q_UNUSED(molecule);

    return QVariant();
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new molecular descriptor.
MolecularDescriptor* MolecularDescriptor::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<MolecularDescriptor>(name);
}

/// Returns a list of available molecular descriptors.
std::vector<std::string> MolecularDescriptor::descriptors()
{
    return PluginManager::instance()->pluginClassNames<MolecularDescriptor>();
}

} // end chemkit namespace
