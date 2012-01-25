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

#include "polymerfile.h"

#include <chemkit/foreach.h>
#include <chemkit/polymer.h>

namespace chemkit {

// === PolymerFilePrivate ================================================== //
class PolymerFilePrivate
{
public:
    std::vector<Polymer *> polymers;
    std::vector<Molecule *> ligands;
};

// === PolymerFile ========================================================= //
/// \class PolymerFile polymerfile.h chemkit/polymerfile.h
/// \ingroup chemkit-io
/// \brief The PolymerFile class contains polymers.
///
/// A list of supported polymer file formats is available at:
/// http://wiki.chemkit.org/Features#Polymer_File_Formats
///
/// \see Polymer

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new polymer file.
PolymerFile::PolymerFile()
    : d(new PolymerFilePrivate)
{
}

/// Creates a new polymer file with \p fileName.
PolymerFile::PolymerFile(const std::string &fileName)
    : GenericFile<PolymerFile, PolymerFileFormat>(fileName),
      d(new PolymerFilePrivate)
{
}

/// Destroys the polymer file object.
PolymerFile::~PolymerFile()
{
    foreach(Polymer *polymer, d->polymers)
        delete polymer;

    foreach(Molecule *ligand, d->ligands)
        delete ligand;

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of polymers in the file.
size_t PolymerFile::size() const
{
    return polymerCount();
}

/// Returns \c true if the file contains no polymers.
bool PolymerFile::isEmpty() const
{
    return size() == 0;
}

// --- File Contents ------------------------------------------------------- //
/// Adds a polymer to the file.
///
/// The ownership of the polymer is passed to the file.
void PolymerFile::addPolymer(Polymer *polymer)
{
    d->polymers.push_back(polymer);
}

/// Removes a polymer from the file and deletes it.
bool PolymerFile::removePolymer(Polymer *polymer)
{
    bool found = takePolymer(polymer);

    if(found){
        delete polymer;
    }

    return found;
}

/// Removes a polymer from the file.
///
/// The ownership of the polymer is passed to the caller.
bool PolymerFile::takePolymer(Polymer *polymer)
{
    std::vector<Polymer *>::iterator location = std::find(d->polymers.begin(), d->polymers.end(), polymer);
    if(location == d->polymers.end()){
        return false;
    }

    d->polymers.erase(location);

    return true;
}

/// Returns the polymer at \p index in the file.
Polymer* PolymerFile::polymer(size_t index) const
{
    return d->polymers[index];
}

/// Returns a list of all the polymers in the file.
std::vector<Polymer *> PolymerFile::polymers() const
{
    return d->polymers;
}

/// Returns the number of polymers in the file.
size_t PolymerFile::polymerCount() const
{
    return d->polymers.size();
}

/// Returns \c true if the file contains \p polymer.
bool PolymerFile::contains(const Polymer *polymer) const
{
    return std::find(d->polymers.begin(), d->polymers.end(), polymer) != d->polymers.end();
}

/// Adds \p ligand to the file.
void PolymerFile::addLigand(Molecule *ligand)
{
    d->ligands.push_back(ligand);
}

/// Removes \p ligand from the file and deletes it.
bool PolymerFile::removeLigand(Molecule *ligand)
{
    bool found = takeLigand(ligand);

    if(found){
        delete ligand;
    }

    return found;
}

/// Removes \p ligand from the file.
///
/// The ownership of \p ligand is passed to the caller.
bool PolymerFile::takeLigand(Molecule *ligand)
{
    std::vector<Molecule *>::iterator location = std::find(d->ligands.begin(),
                                                           d->ligands.end(),
                                                           ligand);
    if(location == d->ligands.end()){
        return false;
    }

    d->ligands.erase(location);

    return true;
}

/// Returns the ligand at \p index in the file.
Molecule* PolymerFile::ligand(size_t index)
{
    assert(index < d->ligands.size());

    return d->ligands[index];
}

/// Returns the ligand in the file with \p name. Returns \c 0 if
/// a ligand with \p name is not found.
Molecule* PolymerFile::ligand(const std::string &name)
{
    foreach(Molecule *ligand, d->ligands){
        if(ligand->name() == name){
            return ligand;
        }
    }

    return 0;
}

/// Returns a list of all the ligands in the file.
std::vector<Molecule *> PolymerFile::ligands() const
{
    return d->ligands;
}

/// Returns the number of ligands in the file.
size_t PolymerFile::ligandCount() const
{
    return d->ligands.size();
}

/// Returns \c true if the file contains \p ligand.
bool PolymerFile::contains(const Molecule *ligand) const
{
    return std::find(d->ligands.begin(),
                     d->ligands.end(),
                     ligand) != d->ligands.end();
}

/// Removes all the polymers and ligands from the file.
void PolymerFile::clear()
{
    foreach(Polymer *polymer, d->polymers)
        delete polymer;

    d->polymers.clear();

    foreach(Molecule *ligand, d->ligands)
        delete ligand;

    d->ligands.clear();
}

} // end chemkit namespace
