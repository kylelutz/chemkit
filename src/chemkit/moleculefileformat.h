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

#ifndef CHEMKIT_MOLECULEFILEFORMAT_H
#define CHEMKIT_MOLECULEFILEFORMAT_H

#include "chemkit.h"

#include <string>
#include <vector>

#include "moleculefile.h"

namespace chemkit {

class MoleculeFileFormatPrivate;

class CHEMKIT_EXPORT MoleculeFileFormat
{
    public:
        // typedefs
        typedef MoleculeFileFormat* (*CreateFunction)();

        // construction and destruction
        virtual ~MoleculeFileFormat();

        // properties
        std::string name() const;

        // options
        void setOption(const std::string &name, const QVariant &value);
        QVariant option(const std::string &name) const;

        // input and output
        virtual bool read(QIODevice *iodev, MoleculeFile *file);
        virtual bool write(const MoleculeFile *file, QIODevice *iodev);

        // error handling
        std::string errorString() const;

        // static methods
        static MoleculeFileFormat* create(const std::string &format);
        static std::vector<std::string> formats();

    protected:
        MoleculeFileFormat(const std::string &name);
        void setErrorString(const std::string &error);

    private:
        MoleculeFileFormatPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEFILEFORMAT_H
