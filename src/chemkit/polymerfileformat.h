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

#ifndef CHEMKIT_POLYMERFILEFORMAT_H
#define CHEMKIT_POLYMERFILEFORMAT_H

#include "chemkit.h"

#include <string>
#include <vector>

#include <QtCore>

namespace chemkit {

class PolymerFile;
class PolymerFileFormatPrivate;

class CHEMKIT_EXPORT PolymerFileFormat
{
    public:
        // typedefs
        typedef PolymerFileFormat* (*CreateFunction)();

        // construction and destruction
        virtual ~PolymerFileFormat();

        // properties
        std::string name() const;

        // input and output
        virtual bool read(QIODevice *iodev, PolymerFile *file);
        virtual bool write(const PolymerFile *file, QIODevice *iodev);

        // error handling
        QString errorString() const;

        // static methods
        static PolymerFileFormat* create(const std::string &name);
        static std::vector<std::string> formats();

    protected:
        PolymerFileFormat(const std::string &name);
        void setErrorString(const QString &errorString);

    private:
        PolymerFileFormatPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_POLYMERFILEFORMAT_H
