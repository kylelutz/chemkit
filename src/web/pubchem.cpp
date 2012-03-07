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

#include "pubchem.h"

#include <QtXml>
#include <QtNetwork>

#include "pubchemquery.h"
#include "downloadthread.h"
#include "pubchemquerythread.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

namespace chemkit {

// === PubChemPrivate ====================================================== //
class PubChemPrivate
{
public:
    QUrl url;
    QString errorString;
};

// === PubChem ============================================================= //
/// \class PubChem pubchem.h chemkit/pubchem.h
/// \ingroup chemkit-web
/// \brief The PubChem class provides access to the %PubChem web
///        API.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new PubChem object.
PubChem::PubChem()
    : d(new PubChemPrivate)
{
    d->url = QUrl("http://pubchem.ncbi.nlm.nih.gov/");
}

/// Destroys the PubChem object.
PubChem::~PubChem()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the url to \p url.
void PubChem::setUrl(const QUrl &url)
{
    d->url = url;
}

/// Returns the url.
QUrl PubChem::url() const
{
    return d->url;
}

// --- Downloads ----------------------------------------------------------- //
/// Downloads and returns the molecule with the compound ID \p id.
/// If an error occurs a null pointer is returned.
boost::shared_ptr<Molecule> PubChem::downloadMolecule(const QString &id) const
{
    QScopedPointer<MoleculeFile> file(downloadFile(id));
    if(!file){
        return boost::shared_ptr<Molecule>();
    }

    return file->molecule();
}

/// Downloads and returns the file with the compound ID \p id. If an
/// error occurs \c 0 is returned.
///
/// The ownership of the returned file is passed to the caller.
MoleculeFile* PubChem::downloadFile(const QString &id) const
{
    QByteArray data = downloadFileData(id, "sdf");
    if(data.isEmpty()){
        return 0;
    }

    std::stringstream buffer(std::string(data.constData(), data.size()));

    MoleculeFile *file = new MoleculeFile;
    file->read(buffer, "sdf");

    return file;
}

/// Downloads and returns the file containing the compounds with ID's
/// in the list \p ids. If an error occurs \c 0 is returned.
///
/// The ownership of the file is passed to the caller.
///
/// For example, to download the file containing %PubChem Compounds
/// 1, 2, 3, 42 and 57:
/// \code
/// QStringList ids;
/// ids << "1" << "2" << "3" << "42" << "57";
/// MoleculeFile *file = pubchem.downloadFile(ids);
/// \endcode
MoleculeFile* PubChem::downloadFile(const QStringList &ids) const
{
    QByteArray data = downloadFileData(ids, "sdf");
    if(data.isEmpty()){
        return 0;
    }

    std::stringstream buffer(std::string(data.constData(), data.size()));

    MoleculeFile *file = new MoleculeFile;
    file->read(buffer, "sdf");

    return file;
}

/// Downloads and returns the file data for the compound with ID
/// \p id. If an error occurs an empty QByteArray is returned.
QByteArray PubChem::downloadFileData(const QString &id, const QString &format) const
{
    Q_UNUSED(format);

    QUrl url(QString("%1summary/summary.cgi?cid=%2&disopt=3DDisplaySDF").arg(d->url.toString())
                                                                        .arg(id));

    return DownloadThread::download(url);
}

/// Downloads and returns the file data for the compounds with ID's
/// in the list \p ids. If an error occurs an empty QByteArray is
/// returned.
QByteArray PubChem::downloadFileData(const QStringList &ids, const QString &format) const
{
    PubChemQuery query = PubChemQuery::downloadQuery(ids, format);

    PubChemQueryThread thread(query);
    thread.start();
    thread.wait();

    // the response contains a URL where the file can be downloaded
    QByteArray response = thread.response();

    QDomDocument document;
    document.setContent(response, false);

    QDomNodeList nodes = document.elementsByTagName("PCT-Download-URL_url");
    if(nodes.isEmpty()){
        return QByteArray();
    }

    QDomNode node = nodes.at(0);
    QDomElement element = node.toElement();
    QString url = element.text();

    return DownloadThread::download(QUrl(url));
}

// --- Search -------------------------------------------------------------- //
/// Searches the %PubChem database for \p query and returns a list
/// of matching compound IDs. The returned list of ids can be passed
/// to PubChem::downloadFile() to download the molecules.
QStringList PubChem::search(const QString &query) const
{
    QString url = "http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term=%1";

    QByteArray response = DownloadThread::download(QUrl(url.arg(query)));

    QDomDocument document;
    document.setContent(response, false);

    QDomNodeList nodes = document.elementsByTagName("Id");

    QStringList ids;

    for(int i = 0; i < nodes.size(); i++){
        QDomNode node = nodes.item(i);
        QDomElement element = node.toElement();
        ids.append(element.text());
    }

    return ids;
}

// --- Standardization ----------------------------------------------------- //
/// Returns a string containing the standardized version of
/// \p formula in \p format. If an error occurs an empty string is
/// returned.
///
/// For example, to standardize a SMILES formula:
/// \code
/// std::string formula = pubchem.standardizeFormula("c3cccc3", "smiles");
/// \endcode
std::string PubChem::standardizeFormula(const std::string &formula, const std::string &format) const
{
    return standardizeFormula(formula, format, format);
}

/// Returns a string containing the standardized version of
/// \p formula from \p inputFormat in \p outputFormat. If an error
/// occurs an empty string is returned.
///
/// For example, to convert an InChI string to standardized SMILES:
/// \code
/// std::string formula = pubchem.standardizeFormula("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi", "smiles");
/// \endcode
std::string PubChem::standardizeFormula(const std::string &formula, const std::string &inputFormat, const std::string &outputFormat) const
{
    if(formula.empty()){
        return std::string();
    }

    PubChemQuery query = PubChemQuery::standardizationQuery(formula, inputFormat, outputFormat);

    PubChemQueryThread thread(query);
    thread.start();
    thread.wait();

    QByteArray response = thread.response();

    QDomDocument document;
    document.setContent(response, false);

    QDomNodeList nodes = document.elementsByTagName("PCT-Structure_structure_string");
    if(nodes.isEmpty()){
        return std::string();
    }

    QDomNode node = nodes.at(0);
    QDomElement element = node.toElement();
    QByteArray elementText = element.text().toAscii();
    std::string standardizedFormula = elementText.constData();

    return standardizedFormula;
}

/// Returns a string containing the standardized formula in \p format
/// for the \p molecule. If an error occurs an empty string is
/// returned.
///
/// For example, to get the standardized InChI formula for a
/// molecule:
/// \code
/// std::string formula = pubchem.standardizeFormula(molecule, "inchi");
/// \endcode
std::string PubChem::standardizeFormula(const Molecule *molecule, const std::string &format) const
{
    return standardizeFormula(molecule->formula("smiles"), "smiles", format);
}

// --- Error Handling ------------------------------------------------------ //
void PubChem::setErrorString(const QString &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occurred.
QString PubChem::errorString() const
{
    return d->errorString;
}

} // end chemkit namespace
