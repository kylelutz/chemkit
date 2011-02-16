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

#ifndef CHEMKIT_BIOCHEMICALFILE_H
#define CHEMKIT_BIOCHEMICALFILE_H

#include "chemkit.h"

#include "protein.h"
#include "nucleicacid.h"
#include "biochemicalfileformat.h"

namespace chemkit {

class BiochemicalFilePrivate;

class CHEMKIT_EXPORT BiochemicalFile
{
    public:
        // construction and destruction
        BiochemicalFile();
        BiochemicalFile(const QString &fileName);
        ~BiochemicalFile();

        // properties
        void setFileName(const QString &fileName);
        QString fileName() const;
        void setFormat(BiochemicalFileFormat *format);
        bool setFormat(const QString &name);
        BiochemicalFileFormat* format();
        const BiochemicalFileFormat* format() const;
        QString formatName() const;

        // file contents
        void addProtein(Protein *protein);
        bool removeProtein(Protein *protein);
        bool deleteProtein(Protein *protein);
        Protein* protein(int index = 0) const;
        QList<Protein *> proteins() const;
        int proteinCount() const;
        void addNucleicAcid(NucleicAcid *nucleicAcid);
        bool removeNucleicAcid(NucleicAcid *nucleicAcid);
        bool deleteNucleicAcid(NucleicAcid *nucleicAcid);
        NucleicAcid* nucleicAcid(int index = 0) const;
        QList<NucleicAcid *> nucleicAcids() const;
        int nucleicAcidCount() const;
        bool contains(const Protein *protein) const;
        bool contains(const NucleicAcid *nucleicAcid) const;

        // input and output
        bool read();
        bool read(const QString &fileName);
        bool read(const QString &fileName, const QString &format);
        bool read(QIODevice *iodev, const QString &format);
        bool write();
        bool write(const QString &fileName);
        bool write(const QString &fileName, const QString &format);
        bool write(QIODevice *iodev);
        bool write(QIODevice *iodev, const QString &format);

        // error handling
        QString errorString() const;

        // static methods
        static QStringList formats();

    private:
        void setErrorString(const QString &error);

    private:
        BiochemicalFilePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_BIOCHEMICALFILE_H
