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

#include "proteindatabank.h"

#include <QtNetwork>

#include "downloadthread.h"

#include <chemkit/polymerfile.h>
#include <chemkit/moleculefile.h>

namespace chemkit {
namespace web {

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
    QScopedPointer<chemkit::io::PolymerFile> file(downloadFile(id));
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

    std::stringstream buffer(std::string(data.constData(), data.size()));

    chemkit::io::MoleculeFile file;
    bool ok = file.read(buffer, "sdf");

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
chemkit::io::PolymerFile* ProteinDataBank::downloadFile(const QString &id) const
{
    QByteArray data = downloadFileData(id, "pdb");
    if(data.isEmpty()){
        return 0;
    }

    std::stringstream buffer(std::string(data.constData(), data.size()));

    chemkit::io::PolymerFile *file = new chemkit::io::PolymerFile;
    file->read(buffer, "pdb");

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
}

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

} // end web namespace
} // end chemkit namespace
