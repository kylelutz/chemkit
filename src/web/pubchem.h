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

#ifndef CHEMKIT_PUBCHEM_H
#define CHEMKIT_PUBCHEM_H

#include "web.h"

#include <QtCore>
#include <QtNetwork>

namespace chemkit {

class Molecule;
class ChemicalFile;
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
        Molecule* downloadMolecule(const QString &id) const;
        ChemicalFile* downloadFile(const QString &id) const;
        ChemicalFile* downloadFile(const QStringList &ids) const;
        QByteArray downloadFileData(const QString &id, const QString &format) const;
        QByteArray downloadFileData(const QStringList &ids, const QString &format) const;

        // search
        QStringList search(const QString &query) const;

        // standardization
        QString standardizeFormula(const QString &formula, const QString &format) const;
        QString standardizeFormula(const QString &formula, const QString &inputFormat, const QString &outputFormat) const;
        QString standardizeFormula(const Molecule *molecule, const QString &format) const;

        // error handling
        QString errorString() const;

    private:
        void setErrorString(const QString &error);

    private:
        PubChemPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PUBCHEM_H
