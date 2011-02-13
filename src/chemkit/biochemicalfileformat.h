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

#ifndef CHEMKIT_BIOCHEMICALFILEFORMAT_H
#define CHEMKIT_BIOCHEMICALFILEFORMAT_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class BiochemicalFile;
class BiochemicalFileFormatPrivate;

class CHEMKIT_EXPORT BiochemicalFileFormat
{
    public:
        // typedefs
        typedef BiochemicalFileFormat* (*CreateFunction)();

        // construction and destruction
        virtual ~BiochemicalFileFormat();

        // properties
        QString name() const;

        // input and output
        virtual bool read(QIODevice *iodev, BiochemicalFile *file);
        virtual bool write(const BiochemicalFile *file, QIODevice *iodev);

        // error handling
        QString errorString() const;

        // static methods
        static BiochemicalFileFormat* create(const QString &format);
        static QStringList formats();
        static void registerFormat(const QString &name, CreateFunction function);
        static void unregisterFormat(const QString &name, CreateFunction function);

    protected:
        BiochemicalFileFormat(const QString &name);
        void setErrorString(const QString &errorString);

    private:
        BiochemicalFileFormatPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_BIOCHEMICALFILEFORMAT_H
