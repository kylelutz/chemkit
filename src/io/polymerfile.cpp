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

#include "moleculefile.h"

#include <chemkit/foreach.h>
#include <chemkit/polymer.h>

namespace chemkit {

// === PolymerFilePrivate ================================================== //
class PolymerFilePrivate
{
public:
    std::vector<boost::shared_ptr<Polymer> > polymers;
    MoleculeFile ligandFile; // use a molecule file to manage ligand molecules
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
/// Adds \p polymer to the file.
void PolymerFile::addPolymer(const boost::shared_ptr<Polymer> &polymer)
{
    d->polymers.push_back(polymer);
}

/// Removes \p polymer from the file.
bool PolymerFile::removePolymer(const boost::shared_ptr<Polymer> &polymer)
{
    std::vector<boost::shared_ptr<Polymer> >::iterator iter = std::find(d->polymers.begin(),
                                                                        d->polymers.end(),
                                                                        polymer);
    if(iter == d->polymers.end()){
        return false;
    }

    d->polymers.erase(iter);
    return true;
}

/// Returns the polymer at \p index in the file.
boost::shared_ptr<Polymer> PolymerFile::polymer(size_t index) const
{
    return d->polymers[index];
}

/// Returns a list of all the polymers in the file.
PolymerFile::PolymerRange PolymerFile::polymers() const
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
    foreach(const boost::shared_ptr<Polymer> &polymerPointer, d->polymers){
        if(polymerPointer.get() == polymer){
            return true;
        }
    }

    return false;
}

/// Adds \p ligand to file.
void PolymerFile::addLigand(const boost::shared_ptr<Molecule> &ligand)
{
    d->ligandFile.addMolecule(ligand);
}

/// Removes \p ligand from the file.
bool PolymerFile::removeLigand(const boost::shared_ptr<Molecule> &ligand)
{
    return d->ligandFile.removeMolecule(ligand);
}

/// Returns the ligand at \p index in the file.
boost::shared_ptr<Molecule> PolymerFile::ligand(size_t index)
{
    return d->ligandFile.molecule(index);
}

/// Returns the ligand in the file with \p name. Returns \c 0 if
/// a ligand with \p name is not found.
boost::shared_ptr<Molecule> PolymerFile::ligand(const std::string &name)
{
    return d->ligandFile.molecule(name);
}

/// Returns a list of all the ligands in the file.
PolymerFile::LigandRange PolymerFile::ligands() const
{
    return d->ligandFile.molecules();
}

/// Returns the number of ligands in the file.
size_t PolymerFile::ligandCount() const
{
    return d->ligandFile.size();
}

/// Returns \c true if the file contains \p ligand.
bool PolymerFile::contains(const boost::shared_ptr<Molecule> &ligand) const
{
    return d->ligandFile.contains(ligand);
}

/// Removes all the polymers and ligands from the file.
void PolymerFile::clear()
{
    d->polymers.clear();
    d->ligandFile.clear();
}

} // end chemkit namespace
