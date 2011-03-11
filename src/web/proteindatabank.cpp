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

#include "proteindatabank.h"

#include <QtNetwork>

#include "downloadthread.h"

#include <chemkit/polymerfile.h>
#include <chemkit/chemicalfile.h>

namespace chemkit {

// === ProteinDataBankPrivate ============================================== //
class ProteinDataBankPrivate
{
	public:
        QUrl url;
        QString errorString;
};

// === ProteinDataBank ===================================================== //
/// \class ProteinDataBank proteindatabank.h chemkit/proteindatabank.h
/// \ingroup chemkit-web
/// \brief The ProteinDataBank class provides access to the %Protein
///        Data Bank.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new protein data bank object.
ProteinDataBank::ProteinDataBank()
	: d(new ProteinDataBankPrivate)
{
    d->url = "http://www.pdb.org/";
}

/// Destroys the protein data bank object.
ProteinDataBank::~ProteinDataBank()
{
	delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the url to \p url.
void ProteinDataBank::setUrl(const QUrl &url)
{
    d->url = url;
}

/// Returns the url.
QUrl ProteinDataBank::url() const
{
    return d->url;
}

// --- Downloads ----------------------------------------------------------- //
/// Downloads and returns the polymer (protein or nucleic acid) with
/// the PDB ID of \p id. If an error occurs \c 0 is returned.
///
/// The ownership of the returned polymer is passed to the caller.
///
/// For example, to download the Lysozyme protein (PDB ID: 2LYZ):
/// \code
/// Polymer *lysozyme = pdb.downloadPolymer("2LYZ");
/// \endcode
Polymer* ProteinDataBank::downloadPolymer(const QString &id) const
{
    QScopedPointer<PolymerFile> file(downloadFile(id));
    if(!file){
        return 0;
    }

    Polymer *polymer = file->polymer();
    file->removePolymer(polymer);

    return polymer;
}

/// Downloads and returns the ligand molecule with \p name. If an
/// error occurs \c 0 is returned.
///
/// The ownership of the returned molecule is passed to the caller.
///
/// For example, to download the heme ligand (named "HEM"):
/// \code
/// Molecule *heme = pdb.downloadLigand("HEM");
/// \endcode
Molecule* ProteinDataBank::downloadLigand(const QString &name) const
{
    QString url("%1/pdb/files/ligand/%2_ideal.sdf");

    QByteArray data = DownloadThread::download(url.arg(d->url.toString()).arg(name.toUpper()));
    if(data.isEmpty()){
        return 0;
    }

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);

    ChemicalFile file;
    bool ok = file.read(&buffer, "sdf");

    if(!ok || file.isEmpty()){
        return 0;
    }

    Molecule *molecule = file.molecule();
    file.removeMolecule(molecule);

    return molecule;
}

/// Downloads the file for the biomolecule with the PDB ID of \p id.
/// If an error occurs \c 0 is returned.
///
/// The ownership of the returned file is passed to the caller.
///
/// For example, to download the ubiquitin pdb file (PDB ID: 1UBQ):
/// \code
/// PolymerFile *file = pdb.downloadFile("1UBQ");
/// \endcode
PolymerFile* ProteinDataBank::downloadFile(const QString &id) const
{
    QByteArray data = downloadFileData(id, "pdb");
    if(data.isEmpty()){
        return 0;
    }

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);

    PolymerFile *file = new PolymerFile;
    file->read(&buffer, "pdb");

    return file;
}

/// Downloads the file data for the biomolecule with the PDB ID of
/// \p id. If an error occurs an empty QByteArray is returned.
QByteArray ProteinDataBank::downloadFileData(const QString &id, const QString &format) const
{
    QUrl url(QString("%1pdb/download/downloadFile.do?fileFormat=%2&compression=NO&structureId=%3").arg(d->url.toString())
                                                                                                  .arg(format.toLower())
                                                                                                  .arg(id.toUpper()));

    return DownloadThread::download(url);
};

// --- Error Handling ------------------------------------------------------ //
void ProteinDataBank::setErrorString(const QString &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
QString ProteinDataBank::errorString() const
{
    return d->errorString;
}

} // end chemkit namespace
