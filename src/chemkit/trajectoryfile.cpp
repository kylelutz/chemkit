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

#include "trajectoryfile.h"

#include "trajectory.h"
#include "trajectoryfileformat.h"

namespace chemkit {

// === TrajectoryFilePrivate =============================================== //
class TrajectoryFilePrivate
{
    public:
        std::string fileName;
        std::string errorString;
        Trajectory *trajectory;
        TrajectoryFileFormat *format;
};

// === TrajectoryFile ====================================================== //
/// \class TrajectoryFile trajectoryfile.h chemkit/trajectoryfile.h
/// \ingroup chemkit
/// \brief The TrajectoryFile class contains a trajectory.
///
/// The following trajectory file formats are supported in chemkit:
///     - \c xtc
///
/// \see Trajectory, TrajectoryFileFormat

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new trajectory file.
TrajectoryFile::TrajectoryFile()
    : d(new TrajectoryFilePrivate)
{
    d->trajectory = 0;
    d->format = 0;
}

/// Creates a new trajectory file with \p fileName.
TrajectoryFile::TrajectoryFile(const std::string &fileName)
    : d(new TrajectoryFilePrivate)
{
    d->fileName = fileName;
    d->trajectory = 0;
    d->format = 0;
}

/// Destroys the trajectory file object.
TrajectoryFile::~TrajectoryFile()
{
    delete d->trajectory;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the file name for the trajectory file to \p fileName.
void TrajectoryFile::setFileName(const std::string &fileName)
{
    d->fileName = fileName;
}

/// Returns the file name for the trajectory file.
std::string TrajectoryFile::fileName() const
{
    return d->fileName;
}

/// Returns \c true if the trajectory file is empty.
bool TrajectoryFile::isEmpty() const
{
    return d->trajectory == 0;
}

// --- File Contents ------------------------------------------------------- //
/// Sets the trajectory for the file to \p trajectory.
void TrajectoryFile::setTrajectory(Trajectory *trajectory)
{
    // delete current trajectory
    deleteTrajectory();

    d->trajectory = trajectory;
}

/// Returns the trajectory that the file contains.
Trajectory* TrajectoryFile::trajectory() const
{
    return d->trajectory;
}

/// Removes the trajectory from the file.
bool TrajectoryFile::removeTrajectory()
{
    if(!d->trajectory){
        return false;
    }

    d->trajectory = 0;
    return true;
}

/// Remove the trajectory from the file and deletes it.
bool TrajectoryFile::deleteTrajectory()
{
    if(!d->trajectory){
        return false;
    }

    delete d->trajectory;
    d->trajectory = 0;
    return true;
}

// --- Input and Output ---------------------------------------------------- //
/// Reads the file.
bool TrajectoryFile::read()
{
    if(d->fileName.empty()){
        return false;
    }

    return read(fileName());
}

/// Reads the file from \p fileName.
bool TrajectoryFile::read(const std::string &fileName)
{
    std::string format = QFileInfo(fileName.c_str()).suffix().toStdString();

    return read(fileName, format);
}

/// Reads the file from \p fileName using format.
bool TrajectoryFile::read(const std::string &fileName, const std::string &format)
{
    QFile file(fileName.c_str());
    if(!file.open(QIODevice::ReadOnly)){
        setErrorString(QString("Failed to open '%1' for reading: %2").arg(fileName.c_str()).arg(file.errorString()).toStdString());
        return false;
    }

    return read(&file, format);
}

/// Reads the file from \p iodev using \p format.
bool TrajectoryFile::read(QIODevice *iodev, const std::string &format)
{
    if(d->format == 0 || d->format->name() != format){
        d->format = TrajectoryFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format.c_str()).toStdString());
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
bool TrajectoryFile::write()
{
    return write(fileName());
}

/// Writes the file to \p fileName.
bool TrajectoryFile::write(const std::string &fileName)
{
    std::string format = QFileInfo(fileName.c_str()).suffix().toStdString();

    return write(fileName, format);
}

/// Writes the file to \p fileName using \p format.
bool TrajectoryFile::write(const std::string &fileName, const std::string &format)
{
    QFile file(fileName.c_str());
    if(!file.open(QIODevice::WriteOnly)){
        setErrorString(QString("Failed to open '%1' for writing: %2").arg(fileName.c_str()).arg(file.errorString()).toStdString());
        return false;
    }

    return write(&file, format);
}

/// Writes the file to \p iodev.
bool TrajectoryFile::write(QIODevice *iodev)
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
bool TrajectoryFile::write(QIODevice *iodev, const std::string &format)
{
    if(!d->format || d->format->name() != format){
        d->format = TrajectoryFileFormat::create(format);
        if(!d->format){
            setErrorString(QString("Format '%1' is not supported").arg(format.c_str()).toStdString());
            iodev->close();
            return false;
        }
    }

    return write(iodev);
}

// --- Error Handling ------------------------------------------------------ //
/// Sets a string describing the last error that occurred.
void TrajectoryFile::setErrorString(const std::string &errorString)
{
    d->errorString = errorString;
}

/// Returns a string describing the last error that occurred.
std::string TrajectoryFile::errorString() const
{
    return d->errorString;
}

// --- Static Methods ------------------------------------------------------ //
/// Returns a list of supported trajectory file formats.
std::vector<std::string> TrajectoryFile::formats()
{
    return TrajectoryFileFormat::formats();
}

} // end chemkit namespace
