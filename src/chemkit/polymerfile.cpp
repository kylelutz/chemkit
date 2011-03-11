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

#include "polymerfile.h"

#include "polymer.h"
#include "polymerfileformat.h"

namespace chemkit {

// === PolymerFilePrivate ================================================== //
class PolymerFilePrivate
{
    public:
        QString fileName;
        QString errorString;
        PolymerFileFormat *format;
        QList<Polymer *> polymers;
};

// === PolymerFile ========================================================= //
/// \class PolymerFile polymerfile.h chemkit/polymerfile.h
/// \ingroup chemkit
/// \brief The PolymerFile class contains polymers.
///
/// \see Polymer

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new polymer file.
PolymerFile::PolymerFile()
    : d(new PolymerFilePrivate)
{
    d->format = 0;
}

/// Creates a new polymer file with \p fileName.
PolymerFile::PolymerFile(const QString &fileName)
    : d(new PolymerFilePrivate)
{
    d->format = 0;
    d->fileName = fileName;
}

/// Destroys the polymer file object.
PolymerFile::~PolymerFile()
{
    qDeleteAll(d->polymers);
    delete d->format;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the file name for the file to \p fileName.
void PolymerFile::setFileName(const QString &fileName)
{
    d->fileName = fileName;
}

/// Returns the file name for the file.
QString PolymerFile::fileName() const
{
    return d->fileName;
}

/// Sets the format for the file to \p format.
void PolymerFile::setFormat(PolymerFileFormat *format)
{
    d->format = format;
}

/// Sets the format for the file to \p name. Returns \c false if an
/// error occured.
bool PolymerFile::setFormat(const QString &name)
{
    PolymerFileFormat *format = PolymerFileFormat::create(name);
    if(!format){
        return false;
    }

    setFormat(format);
    return true;
}

/// Returns the file format for the file.
PolymerFileFormat* PolymerFile::format() const
{
    return d->format;
}

/// Returns the name of the file format for the file.
QString PolymerFile::formatName() const
{
    if(d->format){
        return d->format->name();
    }

    return QString();
}

/// Returns the number of polymers in the file.
int PolymerFile::size() const
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
    d->polymers.append(polymer);
}

/// Removes a polymer from the file.
///
/// The ownership of the polymer is passed to the caller.
bool PolymerFile::removePolymer(Polymer *polymer)
{
    return d->polymers.removeOne(polymer);
}

/// Removes a polymer from the file and deletes it.
bool PolymerFile::deletePolymer(Polymer *polymer)
{
    bool found = removePolymer(polymer);

    if(found){
        delete polymer;
    }

    return found;
}

/// Returns the polymer at \p index in the file.
Polymer* PolymerFile::polymer(int index) const
{
    return d->polymers.value(index, 0);
}

/// Returns a list of all the polymers in the file.
QList<Polymer *> PolymerFile::polymers() const
{
    return d->polymers;
}

/// Returns the number of polymers in the file.
int PolymerFile::polymerCount() const
{
    return d->polymers.size();
}

/// Returns \c true if the file contains \p polymer.
bool PolymerFile::contains(const Polymer *polymer) const
{
    return d->polymers.contains(const_cast<Polymer *>(polymer));
}

/// Removes all the polymers from the file.
void PolymerFile::clear()
{
    qDeleteAll(d->polymers);
    d->polymers.clear();
}

// --- Input and Output ---------------------------------------------------- //
/// Reads the file.
bool PolymerFile::read()
{
    if(d->fileName.isEmpty()){
        return false;
    }

    return read(fileName());
}

/// Reads the file from \p fileName.
bool PolymerFile::read(const QString &fileName)
{
    QString format = QFileInfo(fileName).suffix();

    return read(fileName, format);
}

/// Reads the file from \p fileName using format.
bool PolymerFile::read(const QString &fileName, const QString &format)
{
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly)){
        setErrorString(QString("Failed to open '%1' for reading: %2").arg(fileName).arg(file.errorString()));
        return false;
    }

    return read(&file, format);
}

/// Reads the file from \p iodev using \p format.
bool PolymerFile::read(QIODevice *iodev, const QString &format)
{
    if(d->format == 0 || d->format->name() != format){
        d->format = PolymerFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format));
            iodev->close();
            return false;
        }
    }

    bool ok = d->format->read(iodev, this);
    if(!ok){
        setErrorString(d->format->errorString());
    }

    iodev->close();
    return ok;
}

/// Writes the file.
bool PolymerFile::write()
{
    return write(fileName());
}

/// Writes the file to \p fileName.
bool PolymerFile::write(const QString &fileName)
{
    QString format = QFileInfo(fileName).suffix();

    return write(fileName, format);
}

/// Writes the file to \p fileName using \p format.
bool PolymerFile::write(const QString &fileName, const QString &format)
{
    QFile file(fileName);
    if(!file.open(QIODevice::WriteOnly)){
        setErrorString(QString("Failed to open '%1' for writing: %2").arg(fileName).arg(file.errorString()));
        return false;
    }

    return write(&file, format);
}

/// Writes the file to \p iodev.
bool PolymerFile::write(QIODevice *iodev)
{
    if(!d->format)
        return false;

    bool ok = d->format->write(this, iodev);
    if(!ok){
        setErrorString(d->format->errorString());
    }

    iodev->close();
    return ok;
}

/// Writes the file to \p iodev using \p format.
bool PolymerFile::write(QIODevice *iodev, const QString &format)
{
    if(!d->format || d->format->name() != format){
        d->format = PolymerFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format));
            iodev->close();
            return false;
        }
    }

    return write(iodev);
}

// --- Error Handling ------------------------------------------------------ //
void PolymerFile::setErrorString(const QString &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occured.
QString PolymerFile::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a list of all supported polymer file formats.
QStringList PolymerFile::formats()
{
    return PolymerFileFormat::formats();
}

} // end chemkit namespace
