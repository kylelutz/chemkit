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

#ifndef CHEMKIT_COMMAINITIALIZER_H
#define CHEMKIT_COMMAINITIALIZER_H

#include "chemkit.h"

namespace chemkit {

template<typename T>
class CommaInitializer
{
    public:
        // construction and destruction
        CommaInitializer(T *data, int rowCount, int columnCount);

        // operators
        CommaInitializer<T>& operator,(const T value);

    private:
        T *m_data;
        int m_row;
        int m_column;
        int m_rowCount;
        int m_columnCount;
};

} // end chemkit namespace

#include "commainitializer-inline.h"

#endif // CHEMKIT_COMMAINITIALIZER_H
