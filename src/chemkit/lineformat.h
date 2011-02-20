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
        QString name() const;

        // options
        void setOption(const QString &name, const QVariant &value);
        QVariant option(const QString &name) const;

        // input and output
        virtual bool read(const QString &formula, Molecule *molecule);
        Molecule* read(const QString &formula);
        virtual QString write(const Molecule *molecule);

        // error handling
        QString errorString() const;

        // static methods
        static LineFormat *create(const QString &name);
        static QStringList formats();
        static void registerFormat(const QString &name, CreateFunction function);
        static void unregisterFormat(const QString &name, CreateFunction function);

    protected:
        LineFormat(const QString &name);
        void setErrorString(const QString &error);
        virtual QVariant defaultOption(const QString &name) const;

    private:
        LineFormatPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_LINEFORMAT_H
