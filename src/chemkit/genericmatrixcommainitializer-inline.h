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

#ifndef CHEMKIT_GENERICMATRIXCOMMAINITIALIZER_INLINE_H
#define CHEMKIT_GENERICMATRIXCOMMAINITIALIZER_INLINE_H

#include "genericmatrixcommainitializer.h"

namespace chemkit {

// === GenericMatrixCommaInitializer ======================================= //
// --- Construction and Destruction ---------------------------------------- //
template<typename T>
inline GenericMatrixCommaInitializer<T>::GenericMatrixCommaInitializer(T *data, int rowCount, int columnCount)
    : m_data(data),
      m_rowCount(rowCount),
      m_columnCount(columnCount)
{
    m_row = 0;
    m_column = 0;
}

// --- Operators ----------------------------------------------------------- //
template<typename T>
inline GenericMatrixCommaInitializer<T>& GenericMatrixCommaInitializer<T>::operator,(const T value)
{
    m_column++;
    if(m_column == m_columnCount){
        m_column = 0;
        m_row++;

        // we've filled all the matrix data, so ignore any other values
        if(m_row == m_rowCount){
            return *this;
        }
    }

    // set the value
    m_data[m_column * m_rowCount + m_row] = value;

    return *this;
}

} // end chemkit namespace

#endif // CHEMKIT_GENERICMATRIXCOMMAINITIALIZER_INLINE_H
