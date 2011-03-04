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

#include <chemkit/atomtyper.h>

class SybylAtomTyper : public chemkit::AtomTyper
{
    public:
        // construction and destruction
        SybylAtomTyper();
        ~SybylAtomTyper();

        // types
        QString typeString(int index) const;
        QString typeString(const chemkit::Atom *atom) const;

    protected:
        QString atomType(const chemkit::Atom *atom) const;
        void assignTypes(const chemkit::Molecule *molecule);

    private:
        QVector<QString> m_types;
};
