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

#ifndef CHEMKIT_LINEFORMAT_H
#define CHEMKIT_LINEFORMAT_H

#include "chemkit.h"

#include <string>

#include "molecule.h"

namespace chemkit {

class LineFormatPrivate;

class CHEMKIT_EXPORT LineFormat
{
    public:
        // typedefs
        typedef LineFormat* (*CreateFunction)();

        // construction and destruction
        virtual ~LineFormat();

        // properties
        std::string name() const;

        // options
        void setOption(const std::string &name, const QVariant &value);
        QVariant option(const std::string &name) const;

        // input and output
        virtual bool read(const std::string &formula, Molecule *molecule);
        Molecule* read(const std::string &formula);
        virtual std::string write(const Molecule *molecule);

        // error handling
        std::string errorString() const;

        // static methods
        static LineFormat *create(const std::string &name);
        static QList<std::string> formats();
        static void registerFormat(const std::string &name, CreateFunction function);
        static void unregisterFormat(const std::string &name, CreateFunction function);

    protected:
        LineFormat(const std::string &name);
        void setErrorString(const std::string &error);
        virtual QVariant defaultOption(const std::string &name) const;

    private:
        LineFormatPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_LINEFORMAT_H
