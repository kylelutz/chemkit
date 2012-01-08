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
/// A list of supported molecular descriptors is available at:
/// http://wiki.chemkit.org/Features#Molecular_Descriptors
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
Variant MolecularDescriptor::value(const Molecule *molecule) const
{
    CHEMKIT_UNUSED(molecule);

    return Variant();
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
