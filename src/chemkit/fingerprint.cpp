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

#include "fingerprint.h"

#include "pluginmanager.h"

namespace chemkit {

// === Fingerprint ========================================================= //
/// \class Fingerprint fingerprint.h chemkit/fingerprint.h
/// \ingroup chemkit
/// \brief The Fingerprint class represents a molecular fingerprint.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new fingerprint with \p name.
Fingerprint::Fingerprint(const std::string &name)
    : m_name(name)
{
}

/// Destroys the fingerprint object.
Fingerprint::~Fingerprint()
{
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name for the fingerprint.
std::string Fingerprint::name() const
{
    return m_name;
}

// --- Fingerprint --------------------------------------------------------- //
/// Returns the fingerprint value as a bitset.
Bitset Fingerprint::value(const Molecule *molecule) const
{
    CHEMKIT_UNUSED(molecule);

    return Bitset();
}

// --- Similarity ---------------------------------------------------------- //
/// Returns the tanimoto coefficent between \p a and \p b.
Real Fingerprint::tanimotoCoefficient(const Bitset &a, const Bitset &b)
{
    size_t intersection = (a & b).count();

    return Real(intersection) / Real(a.count() + b.count() - intersection);
}

// --- Static Methods ------------------------------------------------------ //
/// Creates a new fingerprint object for \p name. Returns \c 0 if
/// \p name is not supported.
Fingerprint* Fingerprint::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<Fingerprint>(name);
}

/// Returns a list containing the names of all supported
/// fingerprints.
std::vector<std::string> Fingerprint::fingerprints()
{
    return PluginManager::instance()->pluginClassNames<Fingerprint>();
}

} // end chemkit namespace
