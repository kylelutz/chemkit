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

#include "biochemicalfile.h"

#include "biochemicalfileformat.h"

namespace chemkit {

// === BiochemicalFilePrivate ============================================== //
class BiochemicalFilePrivate
{
    public:
        QString fileName;
        QString errorString;
        BiochemicalFileFormat *format;
        QList<Protein *> proteins;
        QList<NucleicAcid *> nucleicAcids;
};

// === BiochemicalFile ===================================================== //
/// \class BiochemicalFile biochemicalfile.h chemkit/biochemicalfile.h
/// \ingroup chemkit
/// \brief The BiochemicalFile class represents a biochemical data
///        file that contains biomolecules such as proteins and
///        nucleic acids.
///
/// The following biochemical file formats are supported in chemkit:
///     - \c pdb
///     - \c pdbml
///
/// \see ChemicalFile

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty biochemical file object.
BiochemicalFile::BiochemicalFile()
    : d(new BiochemicalFilePrivate)
{
    d->format = 0;
}

/// Creates a new, empty biochemical file object with fileName.
BiochemicalFile::BiochemicalFile(const QString &fileName)
    : d(new BiochemicalFilePrivate)
{
    d->format = 0;
    d->fileName = fileName;
}

/// Destroys the biochemical file object. This also destroys all of
/// the proteins and nucleic acids that it contains.
BiochemicalFile::~BiochemicalFile()
{
    qDeleteAll(d->proteins);
    qDeleteAll(d->nucleicAcids);
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name of the file to \p fileName.
void BiochemicalFile::setFileName(const QString &fileName)
{
    d->fileName = fileName;
}

/// Returns the name of the file.
QString BiochemicalFile::fileName() const
{
    return d->fileName;
}

/// Sets the biochemical file format to \p format.
void BiochemicalFile::setFormat(BiochemicalFileFormat *format)
{
    d->format = format;
}

/// Sets the biochemical file format. Returns \c false if the format
/// \p name is not supported.
bool BiochemicalFile::setFormat(const QString &name)
{
    BiochemicalFileFormat *format = BiochemicalFileFormat::create(name);
    if(!format)
        return false;

    setFormat(format);
    return true;
}

/// Returns the biochemical file format for the file.
BiochemicalFileFormat* BiochemicalFile::format()
{
    return d->format;
}

/// \overload
const BiochemicalFileFormat* BiochemicalFile::format() const
{
    return d->format;
}

/// Returns the name of the biochemical file format.
QString BiochemicalFile::formatName() const
{
    return format()->name();
}

// --- File Contents ------------------------------------------------------- //
/// Adds \p protein to the file.
///
/// The file will take ownership of the protein until it is removed.
void BiochemicalFile::addProtein(Protein *protein)
{
    d->proteins.append(protein);
}

/// Removes \p protein from the file. Returns \c true if \p protein
/// is found and removed successfully.
///
/// The ownership of \p protein is passed to the caller.
bool BiochemicalFile::removeProtein(Protein *protein)
{
    return d->proteins.removeOne(protein);
}

/// Removes \p protein from the file and deletes it. Returns \c true
/// if \p protein is found and deleted successfully.
bool BiochemicalFile::deleteProtein(Protein *protein)
{
    bool found = removeProtein(protein);

    if(found){
        delete protein;
    }

    return found;
}

/// Returns the protein at the \p index.
Protein* BiochemicalFile::protein(int index)
{
    return d->proteins.value(index, 0);
}

/// \overload
const Protein* BiochemicalFile::protein(int index) const
{
    return d->proteins.value(index, 0);
}

/// Returns a list of all the proteins in the file.
QList<Protein *> BiochemicalFile::proteins()
{
    return d->proteins;
}

/// \overload
QList<const Protein *> BiochemicalFile::proteins() const
{
    QList<const Protein *> proteins;

    foreach(const Protein *protein, const_cast<BiochemicalFile *>(this)->proteins()){
        proteins.append(protein);
    }

    return proteins;
}

/// Returns the number of proteins in the file.
int BiochemicalFile::proteinCount() const
{
    return proteins().size();
}

/// Adds \p nucleicAcid to the file.
///
/// The file will take ownership of the nucleic acid until it is
/// removed.
void BiochemicalFile::addNucleicAcid(NucleicAcid *nucleicAcid)
{
    d->nucleicAcids.append(nucleicAcid);
}

/// Removes \p nucleicAcid from the file. Returns \c true if
/// \p nucleicAcid is found and removed successfully.
///
/// The ownership of \p nucleicAcid is passed to the caller.
bool BiochemicalFile::removeNucleicAcid(NucleicAcid *nucleicAcid)
{
    return d->nucleicAcids.removeOne(nucleicAcid);
}

/// Removes \p nucleicAcid from the file and deletes it. Returns
/// \c true if \p nucleicAcid is found and deleted successfully.
bool BiochemicalFile::deleteNucleicAcid(NucleicAcid *nucleicAcid)
{
    bool found = removeNucleicAcid(nucleicAcid);

    if(found){
        delete nucleicAcid;
    }

    return found;
}

/// Returns the nucleic acid at the given index.
NucleicAcid* BiochemicalFile::nucleicAcid(int index)
{
    return d->nucleicAcids.value(index, 0);
}

/// \overload
const NucleicAcid* BiochemicalFile::nucleicAcid(int index) const
{
    return d->nucleicAcids.value(index, 0);
}

/// Returns a list of all the nucleic acids in the file.
QList<NucleicAcid *> BiochemicalFile::nucleicAcids()
{
    return d->nucleicAcids;
}

/// \overload
QList<const NucleicAcid *> BiochemicalFile::nucleicAcids() const
{
    QList<const NucleicAcid *> nucleicAcids;

    foreach(const NucleicAcid *nucleicAcid, const_cast<BiochemicalFile *>(this)->nucleicAcids()){
        nucleicAcids.append(nucleicAcid);
    }

    return nucleicAcids;
}

/// Returns the number of nucleic acids in the file.
int BiochemicalFile::nucleicAcidCount() const
{
    return d->nucleicAcids.size();
}

/// Returns \c true if the file contains the protein.
bool BiochemicalFile::contains(const Protein *protein) const
{
    return d->proteins.contains(const_cast<Protein *>(protein));
}

/// Returns \c true if the file contains the nucleic acid.
bool BiochemicalFile::contains(const NucleicAcid *nucleicAcid) const
{
    return d->nucleicAcids.contains(const_cast<NucleicAcid *>(nucleicAcid));
}

// --- Input and Output ---------------------------------------------------- //
/// Reads the file.
bool BiochemicalFile::read()
{
    if(d->format){
        return read(d->fileName, d->format->name());
    }
    else{
        return read(d->fileName);
    }
}

/// Reads the file from \p fileName.
bool BiochemicalFile::read(const QString &fileName)
{
    QString format = QFileInfo(fileName).suffix();

    return read(fileName, format);
}

/// Reads the file from \p fileName using \p format.
bool BiochemicalFile::read(const QString &fileName, const QString &format)
{
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly)){
        setErrorString(QString("Failed to open '%1' for reading: %2").arg(fileName).arg(file.errorString()));
        return false;
    }

    return read(&file, format);
}

/// Reads the file from \p iodev using \p format.
bool BiochemicalFile::read(QIODevice *iodev, const QString &format)
{
    if(d->format == 0 || d->format->name() != format){
        d->format = BiochemicalFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format));
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

/// Write the file.
bool BiochemicalFile::write()
{
    if(d->format){
        return write(fileName(), d->format->name());
    }
    else{
        return write(fileName());
    }
}

/// Write the file to \p fileName.
bool BiochemicalFile::write(const QString &fileName)
{
    QString format = QFileInfo(fileName).suffix();

    return write(fileName, format);
}

/// Write the file to \p fileName using \p format.
bool BiochemicalFile::write(const QString &fileName, const QString &format)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly)){
        setErrorString(QString("Failed to open '%1' for writing: %2").arg(fileName).arg(file.errorString()));
        return false;
    }

    return write(&file, format);
}

/// Write the file to \p iodev.
bool BiochemicalFile::write(QIODevice *iodev)
{
    if(!d->format)
        return false;

    bool ok = d->format->write(this, iodev);
    if(!ok)
        setErrorString(d->format->errorString());

    iodev->close();
    return ok;
}

/// Write the file to \p iodev using \p format.
bool BiochemicalFile::write(QIODevice *iodev, const QString &format)
{
    if(!d->format || d->format->name() != format){
        d->format = BiochemicalFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format));
            iodev->close();
            return false;
        }
    }

    return write(iodev);
}

// --- Error Handling ------------------------------------------------------ //
void BiochemicalFile::setErrorString(const QString &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
QString BiochemicalFile::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a list of supported biochemical file formats.
QStringList BiochemicalFile::formats()
{
    return BiochemicalFileFormat::formats();
}

} // end chemkit namespace
