/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#ifndef CHEMKIT_VARIANT_H
#define CHEMKIT_VARIANT_H

#include "chemkit.h"

#include <string>

namespace chemkit {

class CHEMKIT_EXPORT Variant
{
    public:
        // enumerations
        enum Type {
            Null,
            Bool,
            Int,
            Long,
            Float,
            Double,
            Pointer,
            String
        };

        // construction and destruction
        inline Variant();
        template<typename T> Variant(T value);
        inline Variant(const Variant &variant);
        inline ~Variant();

        // properties
        inline Type type() const;
        inline bool isNull() const;

        // value
        template<typename T> bool setValue(T value);
        template<typename T> T value() const;
        inline void clear();

        // conversions
        inline bool toBool() const;
        inline char toChar() const;
        inline unsigned char toUChar() const;
        inline short toShort() const;
        inline unsigned short toUShort() const;
        inline int toInt() const;
        inline unsigned int toUInt() const;
        inline long toLong() const;
        inline unsigned long toULong() const;
        inline float toFloat() const;
        inline double toDouble() const;
        inline void* toPointer() const;
        inline std::string toString() const;

        // operators
        inline Variant& operator=(const Variant &variant);

    private:
        Type m_type;
        union {
            bool _bool;
            char _char;
            int _int;
            long _long;
            float _float;
            double _double;
            void *pointer;
            std::string *string;
        } m_value;
};

} // end chemkit namespace

#include "variant-inline.h"

#endif // CHEMKIT_VARIANT_H
