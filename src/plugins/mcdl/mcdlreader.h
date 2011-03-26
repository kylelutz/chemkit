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

#ifndef MCDLREADER_H
#define MCDLREADER_H

#include <QtCore>

#include <chemkit/molecule.h>

class McdlReader
{
    public:
        // construction and destruction
        McdlReader();
        ~McdlReader();

        // reading
        bool read(const std::string &formula, chemkit::Molecule *molecule);
        bool read(const char *formula, chemkit::Molecule *molecule);

        // error handling
        std::string errorString() const;

    private:
        bool readCompositionModule();
        bool readConnectionModule();
        int readNumber(const char **p);
        int readElement(const char **p);
        void addFragmentCopies(chemkit::Atom *atom, int quantity);
        void addFragmentConnections(const QList<int> &connections, int fragment);
        void setErrorString(const std::string &error);

    private:
        const char *p;
        const char *m_formula;
        chemkit::Molecule *m_molecule;
        QList<chemkit::Atom *> m_fragments;
        std::string m_errorString;
};

#endif // MCDLREADER_H
