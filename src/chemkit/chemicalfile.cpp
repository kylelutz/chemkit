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

#include "chemicalfile.h"

#include "molecule.h"
#include "chemicalfileformat.h"

namespace chemkit {

// === ChemicalFilePrivate ================================================= //
class ChemicalFilePrivate
{
    public:
        std::string fileName;
        std::string errorString;
        ChemicalFileFormat *format;
        QList<Molecule *> molecules;
        std::map<std::string, QVariant> fileData;
        std::map<const Molecule *, std::map<std::string, QVariant> > moleculeData;
};

// === ChemicalFile ======================================================== //
/// \class ChemicalFile chemicalfile.h chemkit/chemicalfile.h
/// \ingroup chemkit
/// \brief The ChemicalFile class represents a chemical data file
///        containing one or more molecules.
///
/// Chemical files object can be used to both read and write molecule
/// data contained in files.
///
/// The following chemical file formats are supported in chemkit:
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
/// The following example shows how to read a molecule from a chemical
/// file:
/// \code
/// // create file
/// ChemicalFile file("ethanol.mol");
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
ChemicalFile::ChemicalFile()
    : d(new ChemicalFilePrivate)
{
    d->format = 0;
}

/// Creates a new, empty file object with \p fileName.
ChemicalFile::ChemicalFile(const std::string &fileName)
    : d(new ChemicalFilePrivate)
{
    d->format = 0;
    d->fileName = fileName;
}

/// Destroys the file object. Destroying the file will also destroy
/// any molecules that it contains.
ChemicalFile::~ChemicalFile()
{
    qDeleteAll(d->molecules);
    delete d->format;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name of the file to \p fileName.
void ChemicalFile::setFileName(const std::string &fileName)
{
    d->fileName = fileName;
}

/// Returns the name of the file.
std::string ChemicalFile::fileName() const
{
    return d->fileName;
}

/// Sets the format for the file to \p format.
void ChemicalFile::setFormat(ChemicalFileFormat *format)
{
    d->format = format;
}

/// Sets the format of the file to \p name. If name is not a valid
/// format the current format will remain unchanged and \c false
/// will be returned.
bool ChemicalFile::setFormat(const std::string &name)
{
    ChemicalFileFormat *format = ChemicalFileFormat::create(name);

    if(!format){
        return false;
    }
    else{
        setFormat(format);
        return true;
    }
}

/// Returns the format object for the file.
ChemicalFileFormat* ChemicalFile::format() const
{
    return d->format;
}

/// Returns the name of the format for this file.
std::string ChemicalFile::formatName() const
{
    if(d->format){
        return d->format->name();
    }

    return std::string();
}

/// Returns the number of molecules in the file.
int ChemicalFile::size() const
{
    return moleculeCount();
}

/// Returns \c true if the file contains no molecules (i.e.
/// size() \c == \c 0).
bool ChemicalFile::isEmpty() const
{
    return size() == 0;
}

// --- File Contents ------------------------------------------------------- //
/// Adds the molecule to the file.
///
/// The file will take ownership of the molecule until it is removed.
void ChemicalFile::addMolecule(Molecule *molecule)
{
    d->molecules.append(molecule);
}

/// Removes the molecule from the file. Returns \c true if
/// \p molecule is found and removed successfully.
///
/// The ownership of \p molecule is passed to the caller.
bool ChemicalFile::removeMolecule(Molecule *molecule)
{
    d->moleculeData.erase(molecule);
    return d->molecules.removeOne(molecule);
}

/// Removes the molecule from the file and deletes it. Returns
/// \c true if \p molecule is found and deleted successfully.
bool ChemicalFile::deleteMolecule(Molecule *molecule)
{
    bool found = removeMolecule(molecule);

    if(found){
        delete molecule;
    }

    return found;
}

/// Returns a list of all the molecules in the file.
QList<Molecule *> ChemicalFile::molecules() const
{
    return d->molecules;
}

/// Returns the number of molecules in the file.
int ChemicalFile::moleculeCount() const
{
    return d->molecules.size();
}

/// Returns the molecule at \p index in the file.
Molecule* ChemicalFile::molecule(int index) const
{
    return d->molecules.value(index, 0);
}

/// Returns \c true if the file contains \p molecule.
bool ChemicalFile::contains(const Molecule *molecule) const
{
    return d->molecules.contains(const_cast<Molecule *>(molecule));
}

/// Removes all of the molecules from the file and deletes all
/// of the data in the file.
void ChemicalFile::clear()
{
    qDeleteAll(d->molecules);
    d->molecules.clear();
    d->fileData.clear();
    d->moleculeData.clear();
}

// --- File Data ----------------------------------------------------------- //
/// Sets data with \p name to \p value for the file.
void ChemicalFile::setFileData(const std::string &name, const QVariant &value)
{
    d->fileData[name] = value;
}

/// Returns the data for \p name.
QVariant ChemicalFile::fileData(const std::string &name) const
{
    std::map<std::string, QVariant>::iterator element = d->fileData.find(name);
    if(element != d->fileData.end()){
        return element->second;
    }

    return QVariant();
}

/// Sets data for \p molecule with \p name to \p value in the file.
void ChemicalFile::setMoleculeData(const Molecule *molecule, const std::string &name, const QVariant &value)
{
    d->moleculeData[molecule][name] = value;
}

/// Returns data for \p molecule with \p name in the file.
QVariant ChemicalFile::moleculeData(const Molecule *molecule, const std::string &name) const
{
    const std::map<std::string, QVariant> &moleculeDataMap = d->moleculeData[molecule];
    std::map<std::string, QVariant>::const_iterator element = moleculeDataMap.find(name);
    if(element != moleculeDataMap.end()){
        return element->second;
    }

    return QVariant();
}

// --- Input and Output ---------------------------------------------------- //
/// Reads the file.
bool ChemicalFile::read()
{
    if(d->fileName.empty()){
        return false;
    }

    return read(fileName());
}

/// Reads the file from \p fileName.
bool ChemicalFile::read(const std::string &fileName)
{
    std::string format = QFileInfo(fileName.c_str()).suffix().toStdString();

    return read(fileName, format);
}

/// Reads the file from \p fileName using format.
bool ChemicalFile::read(const std::string &fileName, const std::string &format)
{
    QFile file(fileName.c_str());
    if(!file.open(QIODevice::ReadOnly)){
        setErrorString(QString("Failed to open '%1' for reading: %2").arg(fileName.c_str()).arg(file.errorString()).toStdString());
        return false;
    }

    return read(&file, format);
}

/// Reads the file from \p iodev using \p format.
bool ChemicalFile::read(QIODevice *iodev, const std::string &format)
{
    if(d->format == 0 || d->format->name() != format){
        d->format = ChemicalFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format.c_str()).toStdString());
            iodev->close();
            return false;
        }
    }

    bool ok = d->format->read(iodev, this);
    if(!ok)
        setErrorString(d->format->errorString());

    iodev->close();
    return ok;
}

/// Writes the file.
bool ChemicalFile::write()
{
    return write(fileName());
}

/// Writes the file to \p fileName.
bool ChemicalFile::write(const std::string &fileName)
{
    std::string format = QFileInfo(fileName.c_str()).suffix().toStdString();

    return write(fileName, format);
}

/// Writes the file to \p fileName using \p format.
bool ChemicalFile::write(const std::string &fileName, const std::string &format)
{
    QFile file(fileName.c_str());
    if(!file.open(QIODevice::WriteOnly)){
        setErrorString(QString("Failed to open '%1' for writing: %2").arg(fileName.c_str()).arg(file.errorString()).toStdString());
        return false;
    }

    return write(&file, format);
}

/// Writes the file to \p iodev.
bool ChemicalFile::write(QIODevice *iodev)
{
    if(!d->format)
        return false;

    bool ok = d->format->write(this, iodev);
    if(!ok)
        setErrorString(d->format->errorString());

    iodev->close();
    return ok;
}

/// Writes the file to \p iodev using \p format.
bool ChemicalFile::write(QIODevice *iodev, const std::string &format)
{
    if(!d->format || d->format->name() != format){
        d->format = ChemicalFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format.c_str()).toStdString());
            iodev->close();
            return false;
        }
    }

    return write(iodev);
}

// --- Error Handling ------------------------------------------------------ //
void ChemicalFile::setErrorString(const std::string &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
std::string ChemicalFile::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a list of all supported chemical file formats.
std::vector<std::string> ChemicalFile::formats()
{
    return ChemicalFileFormat::formats();
}

/// Reads and returns a molecule from the file. Returns \c 0 if there
/// was an error reading the file or the file is empty.
///
/// This static convenience method allows for the reading of molecule
/// from a file without explicitly creating a file object.
Molecule* ChemicalFile::quickRead(const std::string &fileName)
{
    ChemicalFile file(fileName);

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
void ChemicalFile::quickWrite(const Molecule *molecule, const std::string &fileName)
{
    ChemicalFile file;
    file.addMolecule(const_cast<Molecule *>(molecule));
    file.write(fileName);

    // remove the molecule from the file so that it does
    // not get deleted when this function returns
    file.removeMolecule(const_cast<Molecule *>(molecule));
}

} // end chemkit namespace
