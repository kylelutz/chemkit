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

#ifndef CHEMKIT_POLYMERFILE_H
#define CHEMKIT_POLYMERFILE_H

#include "chemkit.h"

#include <string>
#include <vector>

#include <QtCore>

namespace chemkit {

class Polymer;
class PolymerFileFormat;
class PolymerFilePrivate;

class CHEMKIT_EXPORT PolymerFile
{
    public:
        // construction and destruction
        PolymerFile();
        PolymerFile(const QString &fileName);
        ~PolymerFile();

        // properties
        void setFileName(const QString &fileName);
        QString fileName() const;
        void setFormat(PolymerFileFormat *format);
        bool setFormat(const std::string &name);
        PolymerFileFormat* format() const;
        std::string formatName() const;
        int size() const;
        bool isEmpty() const;

        // file contents
        void addPolymer(Polymer *polymer);
        bool removePolymer(Polymer *polymer);
        bool deletePolymer(Polymer *polymer);
        Polymer* polymer(int index = 0) const;
        QList<Polymer *> polymers() const;
        int polymerCount() const;
        bool contains(const Polymer *polymer) const;
        void clear();

        // input and output
        bool read();
        bool read(const QString &fileName);
        bool read(const QString &fileName, const std::string &format);
        bool read(QIODevice *iodev, const std::string &format);
        bool write();
        bool write(const QString &fileName);
        bool write(const QString &fileName, const std::string &format);
        bool write(QIODevice *iodev);
        bool write(QIODevice *iodev, const std::string &format);

        // error handling
        QString errorString() const;

        // static methods
        static std::vector<std::string> formats();

    private:
        void setErrorString(const QString &errorString);

    private:
        PolymerFilePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_POLYMERFILE_H
