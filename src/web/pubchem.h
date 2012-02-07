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

#ifndef CHEMKIT_PUBCHEM_H
#define CHEMKIT_PUBCHEM_H

#include "web.h"

#include <string>

#include <boost/shared_ptr.hpp>

#include <QtCore>
#include <QtNetwork>

namespace chemkit {

class Molecule;
class MoleculeFile;
class PubChemPrivate;

class CHEMKIT_WEB_EXPORT PubChem
{
public:
    // construction and destruction
    PubChem();
    ~PubChem();

    // properties
    void setUrl(const QUrl &url);
    QUrl url() const;

    // downloads
    boost::shared_ptr<Molecule> downloadMolecule(const QString &id) const;
    MoleculeFile* downloadFile(const QString &id) const;
    MoleculeFile* downloadFile(const QStringList &ids) const;
    QByteArray downloadFileData(const QString &id, const QString &format) const;
    QByteArray downloadFileData(const QStringList &ids, const QString &format) const;

    // search
    QStringList search(const QString &query) const;

    // standardization
    std::string standardizeFormula(const std::string &formula, const std::string &format) const;
    std::string standardizeFormula(const std::string &formula, const std::string &inputFormat, const std::string &outputFormat) const;
    std::string standardizeFormula(const Molecule *molecule, const std::string &format) const;

    // error handling
    QString errorString() const;

private:
    void setErrorString(const QString &error);

private:
    PubChemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PUBCHEM_H
