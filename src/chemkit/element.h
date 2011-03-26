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

#ifndef CHEMKIT_ELEMENT_H
#define CHEMKIT_ELEMENT_H

#include "chemkit.h"

#include <string>

namespace chemkit {

class CHEMKIT_EXPORT Element
{
    public:
        // construction and destruction
        Element();
        Element(int atomicNumber);
        Element(const char *symbol);
        Element(const std::string &symbol);

        // properties
        void setAtomicNumber(int atomicNumber);
        int atomicNumber() const;
        std::string symbol() const;
        std::string name() const;
        int period() const;
        Float mass() const;
        Float electronegativity() const;
        Float covalentRadius() const;
        Float vanDerWaalsRadius() const;
        int expectedValence() const;
        bool isValid() const;
        bool isMetal() const;
        bool isNonmetal() const;

        // operators
        bool operator==(const Element &element) const;

        // static methods
        static int atomicNumber(const std::string &symbol);
        static int atomicNumber(const char *symbol);
        static int atomicNumber(const char *symbol, int length);
        static int atomicNumber(char symbol);
        static bool isValidAtomicNumber(int atomicNumber);
        static bool isValidSymbol(const std::string &symbol);

    private:
        int m_atomicNumber;
};

} // end chemkit namespace

#include "element-inline.h"

#endif // CHEMKIT_ELEMENT_H
