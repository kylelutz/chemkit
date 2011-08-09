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

#include "moleculefile.h"

#include <chemkit/foreach.h>
#include <chemkit/molecule.h>

namespace chemkit {

// === MoleculeFilePrivate ================================================= //
class MoleculeFilePrivate
{
    public:
        std::vector<Molecule *> molecules;
        std::map<std::string, Variant> fileData;
};

// === MoleculeFile ======================================================== //
/// \class MoleculeFile moleculefile.h chemkit/moleculefile.h
/// \ingroup chemkit-io
/// \brief The MoleculeFile class represents a molecular data file
///        containing one or more molecules.
///
/// Molecule files object can be used to both read and write molecule
/// data contained in files.
///
/// The following molecule file formats are supported in chemkit:
///     - \c cml
///     - \c inchi
///     - \c mdl
///     - \c mol
///     - \c mol2
///     - \c sd
///     - \c sdf
///     - \c smi
///     - \c txyz
///     - \c xyz
///
/// The following example shows how to read a molecule from a file:
/// \code
/// // create file
/// MoleculeFile file("ethanol.mol");
///
/// // read file
/// file.read();
///
/// // access molecule
/// Molecule *molecule = file.molecule();
/// \endcode
///
/// \see PolymerFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty file object.
MoleculeFile::MoleculeFile()
    : d(new MoleculeFilePrivate)
{
}

/// Creates a new, empty file object with \p fileName.
MoleculeFile::MoleculeFile(const std::string &fileName)
    : GenericFile<MoleculeFile, MoleculeFileFormat>(fileName),
      d(new MoleculeFilePrivate)
{
}

/// Destroys the file object. Destroying the file will also destroy
/// any molecules that it contains.
MoleculeFile::~MoleculeFile()
{
    foreach(Molecule *molecule, d->molecules)
        delete molecule;

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of molecules in the file.
int MoleculeFile::size() const
{
    return moleculeCount();
}

/// Returns \c true if the file contains no molecules (i.e.
/// size() \c == \c 0).
bool MoleculeFile::isEmpty() const
{
    return size() == 0;
}

// --- File Contents ------------------------------------------------------- //
/// Adds the molecule to the file.
///
/// The file will take ownership of the molecule until it is removed.
void MoleculeFile::addMolecule(Molecule *molecule)
{
    d->molecules.push_back(molecule);
}

/// Removes the molecule from the file. Returns \c true if
/// \p molecule is found and removed successfully.
///
/// The ownership of \p molecule is passed to the caller.
bool MoleculeFile::removeMolecule(Molecule *molecule)
{
    std::vector<Molecule *>::iterator location = std::find(d->molecules.begin(), d->molecules.end(), molecule);
    if(location == d->molecules.end()){
        return false;
    }

    d->molecules.erase(location);

    return true;
}

/// Removes the molecule from the file and deletes it. Returns
/// \c true if \p molecule is found and deleted successfully.
bool MoleculeFile::deleteMolecule(Molecule *molecule)
{
    bool found = removeMolecule(molecule);

    if(found){
        delete molecule;
    }

    return found;
}

/// Returns a list of all the molecules in the file.
std::vector<Molecule *> MoleculeFile::molecules() const
{
    return d->molecules;
}

/// Returns the number of molecules in the file.
int MoleculeFile::moleculeCount() const
{
    return d->molecules.size();
}

/// Returns the molecule at \p index in the file.
Molecule* MoleculeFile::molecule(int index) const
{
    return d->molecules[index];
}

/// Returns \c true if the file contains \p molecule.
bool MoleculeFile::contains(const Molecule *molecule) const
{
    return std::find(d->molecules.begin(), d->molecules.end(), molecule) != d->molecules.end();
}

/// Removes all of the molecules from the file and deletes all
/// of the data in the file.
void MoleculeFile::clear()
{
    foreach(Molecule *molecule, d->molecules)
        delete molecule;

    d->molecules.clear();
    d->fileData.clear();
}

// --- Static Methods ------------------------------------------------------ //
/// Reads and returns a molecule from the file. Returns \c 0 if there
/// was an error reading the file or the file is empty.
///
/// This static convenience method allows for the reading of molecule
/// from a file without explicitly creating a file object.
Molecule* MoleculeFile::quickRead(const std::string &fileName)
{
    MoleculeFile file(fileName);

    if(!file.read() || file.isEmpty()){
        return 0;
    }

    Molecule *molecule = file.molecule();

    // remove the molecule from the file so that it does
    // not get deleted when this function returns
    file.removeMolecule(molecule);

    return molecule;
}

/// Writes \p molecule to the file with \p fileName.
///
/// This static convenience method allows for the writing of molecule
/// to a file without explicitly creating a file object.
void MoleculeFile::quickWrite(const Molecule *molecule, const std::string &fileName)
{
    MoleculeFile file;
    file.addMolecule(const_cast<Molecule *>(molecule));
    file.write(fileName);

    // remove the molecule from the file so that it does
    // not get deleted when this function returns
    file.removeMolecule(const_cast<Molecule *>(molecule));
}

} // end chemkit namespace
