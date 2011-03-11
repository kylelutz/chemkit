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

#ifndef CHEMKIT_PROTEINDATABANK_H
#define CHEMKIT_PROTEINDATABANK_H

#include "web.h"

#include <QtCore>

namespace chemkit {

class Polymer;
class Molecule;
class PolymerFile;
class ProteinDataBankPrivate;

class CHEMKIT_WEB_EXPORT ProteinDataBank
{
    public:
        // construction and destruction
        ProteinDataBank();
        ~ProteinDataBank();

        // properties
        void setUrl(const QUrl &url);
        QUrl url() const;

        // downloads
        Polymer* downloadPolymer(const QString &id) const;
        Molecule* downloadLigand(const QString &name) const;
        PolymerFile* downloadFile(const QString &id) const;
        QByteArray downloadFileData(const QString &id, const QString &format) const;

        // error handling
        QString errorString() const;

    private:
        void setErrorString(const QString &error);

    private:
        ProteinDataBankPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PROTEINDATABANK_H
