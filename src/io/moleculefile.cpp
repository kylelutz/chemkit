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

#include <boost/lambda/lambda.hpp>

#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/variantmap.h>

namespace chemkit {

// === MoleculeFilePrivate ================================================= //
class MoleculeFilePrivate
{
public:
    std::vector<boost::shared_ptr<Molecule> > molecules;
    VariantMap fileData;
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
/// The MoleculeFile class stores molecules using boost::shared_ptr's
/// which allow the each molecule's memory to be shared between
/// multiple classes.
///
/// A list of supported molecule file formats is available at:
/// http://wiki.chemkit.org/Features#Molecule_File_Formats
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
/// boost::shared_ptr<Molecule> molecule = file.molecule();
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
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the number of molecules in the file.
size_t MoleculeFile::size() const
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
void MoleculeFile::addMolecule(const boost::shared_ptr<Molecule> &molecule)
{
    d->molecules.push_back(molecule);
}

/// Removes the molecule from the file. Returns \c true if
/// \p molecule is found and removed successfully.
bool MoleculeFile::removeMolecule(const boost::shared_ptr<Molecule> &molecule)
{
    std::vector<boost::shared_ptr<Molecule> >::iterator iter = std::find(d->molecules.begin(),
                                                                         d->molecules.end(),
                                                                         molecule);
    if(iter == d->molecules.end()){
        return false;
    }

    d->molecules.erase(iter);
    return true;
}

/// Returns a range containing all of the molecules in the file.
MoleculeFile::MoleculeRange MoleculeFile::molecules() const
{
    return boost::make_iterator_range(d->molecules.begin(),
                                      d->molecules.end());
}

/// Returns the number of molecules in the file.
size_t MoleculeFile::moleculeCount() const
{
    return d->molecules.size();
}

/// Returns the molecule at \p index in the file.
boost::shared_ptr<Molecule> MoleculeFile::molecule(size_t index) const
{
    return d->molecules[index];
}

/// Returns the molecule in the file with \p name. Returns a null
/// pointer if no molecule with \p name is found.
boost::shared_ptr<Molecule> MoleculeFile::molecule(const std::string &name) const
{
    foreach(const boost::shared_ptr<Molecule> &molecule, d->molecules){
        if(molecule->name() == name){
            return molecule;
        }
    }

    return boost::shared_ptr<Molecule>();
}

/// Returns \c true if the file contains \p molecule.
bool MoleculeFile::contains(const boost::shared_ptr<Molecule> &molecule) const
{
    return std::find(d->molecules.begin(),
                     d->molecules.end(),
                     molecule) != d->molecules.end();
}

/// Removes all of the molecules from the file and deletes all
/// of the data in the file.
void MoleculeFile::clear()
{
    d->molecules.clear();
    d->fileData.clear();
}

// --- Static Methods ------------------------------------------------------ //
/// Reads and returns a molecule from the file. Returns a null pointer if
/// there was an error reading the file or the file is empty.
///
/// This static convenience method allows for the reading of molecule
/// from a file without explicitly creating a file object.
boost::shared_ptr<Molecule> MoleculeFile::quickRead(const std::string &fileName)
{
    MoleculeFile file(fileName);

    if(!file.read() || file.isEmpty()){
        return boost::shared_ptr<Molecule>();
    }

    return file.molecule();
}

/// Writes \p molecule to the file with \p fileName.
///
/// This static convenience method allows for the writing of molecule
/// to a file without explicitly creating a file object.
void MoleculeFile::quickWrite(const Molecule *molecule, const std::string &fileName)
{
    // create a shared_ptr with a null deleter
    boost::shared_ptr<Molecule> moleculePointer(const_cast<Molecule *>(molecule),
                                                boost::lambda::_1);

    MoleculeFile file;
    file.addMolecule(moleculePointer);
    file.write(fileName);
}

} // end chemkit namespace
