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

#include "pubchem.h"

#include <QtXml>
#include <QtNetwork>

#include "pubchemquery.h"
#include "downloadthread.h"
#include "pubchemquerythread.h"

#include <chemkit/molecule.h>
#include <chemkit/chemicalfile.h>

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
/// If an error occurs \c 0 is returned.
///
/// The ownership of the returned molecule is passed to the caller.
Molecule* PubChem::downloadMolecule(const QString &id) const
{
    QScopedPointer<ChemicalFile> file(downloadFile(id));
    if(!file){
        return 0;
    }

    Molecule *molecule = file->molecule();
    file->removeMolecule(molecule);

    return molecule;
}

/// Downloads and returns the file with the compound ID \p id. If an
/// error occurs \c 0 is returned.
///
/// The ownership of the returned file is passed to the caller.
ChemicalFile* PubChem::downloadFile(const QString &id) const
{
    QByteArray data = downloadFileData(id, "sdf");
    if(data.isEmpty()){
        return 0;
    }

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);

    ChemicalFile *file = new ChemicalFile;
    file->read(&buffer, "sdf");

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
/// ChemicalFile *file = pubchem.downloadFile(ids);
/// \endcode
ChemicalFile* PubChem::downloadFile(const QStringList &ids) const
{
    QByteArray data = downloadFileData(ids, "sdf");
    if(data.isEmpty()){
        return 0;
    }

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);

    ChemicalFile *file = new ChemicalFile;
    file->read(&buffer, "sdf");

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
/// \p formula in \p format. If an error occurs an empty QString is
/// returned.
///
/// For example, to standardize a SMILES formula:
/// \code
/// QString formula = pubchem.standardizeFormula("c3cccc3", "smiles");
/// \endcode
QString PubChem::standardizeFormula(const QString &formula, const QString &format) const
{
    return standardizeFormula(formula, format, format);
}

/// Returns a string containing the standardized version of
/// \p formula from \p inputFormat in \p outputFormat. If an error
/// occurs an empty QString is returned.
///
/// For example, to convert an InChI string to standardized SMILES:
/// \code
/// QString formula = pubchem.standardizeFormula("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi", "smiles");
/// \endcode
QString PubChem::standardizeFormula(const QString &formula, const QString &inputFormat, const QString &outputFormat) const
{
    if(formula.isEmpty()){
        return QString();
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
        return QString();
    }

    QDomNode node = nodes.at(0);
    QDomElement element = node.toElement();
    QString standardizedFormula = element.text();

    return standardizedFormula;
}

/// Returns a string containing the standardized formula in \p format
/// for the \p molecule. If an error occurs an empty QString is
/// returned.
///
/// For example, to get the standardized InChI formula for a
/// molecule:
/// \code
/// QString formula = pubchem.standardizeFormula(molecule, "inchi");
/// \endcode
QString PubChem::standardizeFormula(const Molecule *molecule, const QString &format) const
{
    return standardizeFormula(molecule->formula("smiles").c_str(), "smiles", format);
}

// --- Error Handling ------------------------------------------------------ //
void PubChem::setErrorString(const QString &error)
{
    d->errorString = error;
}

/// Returns a string describing the last error that occured.
QString PubChem::errorString() const
{
    return d->errorString;
}

} // end chemkit namespace
