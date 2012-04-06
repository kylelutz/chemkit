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

#ifndef CHEMKIT_VARIANT_INLINE_H
#define CHEMKIT_VARIANT_INLINE_H

#include "variant.h"

#include <boost/lexical_cast.hpp>

namespace chemkit {

// === Variant ============================================================= //
/// \class Variant variant.h chemkit/variant.h
/// \ingroup chemkit
/// \brief The Variant class represents a union of data values.
///
/// Variant objects allow for the storage of and conversion between
/// a variety of different data types.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a null variant.
inline Variant::Variant()
    : m_type(Null)
{
}

/// Creates a variant to store \p value.
template<typename T>
inline Variant::Variant(T value)
    : m_type(Null)
{
    setValue(value);
}

/// Creates a new copy of \p variant.
inline Variant::Variant(const Variant &variant)
    : m_type(variant.type())
{
    if(m_type == String){
        m_value.string = new std::string(variant.toString());
    }
    else if(m_type != Null){
        m_value = variant.m_value;
    }
}

/// Destroys the variant object
inline Variant::~Variant()
{
    clear();
}

// --- Properties ---------------------------------------------------------- //
/// Returns variant's type.
inline Variant::Type Variant::type() const
{
    return m_type;
}

/// Returns \c true if the variant is null.
inline bool Variant::isNull() const
{
    return m_type == Null;
}

// --- Value --------------------------------------------------------------- //
/// Sets the value of the variant to \p value.
template<typename T>
inline bool Variant::setValue(T value)
{
    CHEMKIT_UNUSED(value);

    clear();

    return false;
}

template<>
inline bool Variant::setValue(bool value)
{
    clear();

    m_type = Bool;
    m_value._bool = value;

    return true;
}

template<>
inline bool Variant::setValue(char value)
{
    clear();

    m_type = Int;
    m_value._int = value;

    return true;
}

template<>
inline bool Variant::setValue(short value)
{
    clear();

    m_type = Int;
    m_value._int = value;

    return true;
}

template<>
inline bool Variant::setValue(int value)
{
    clear();

    m_type = Int;
    m_value._int = value;

    return true;
}

template<>
inline bool Variant::setValue(long value)
{
    clear();

    m_type = Long;
    m_value._long = value;

    return true;
}

template<>
inline bool Variant::setValue(size_t value)
{
    clear();

    m_type = Long;
    m_value._long = value;

    return true;
}

template<>
inline bool Variant::setValue(float value)
{
    clear();

    m_type = Float;
    m_value._float = value;

    return true;
}

template<>
inline bool Variant::setValue(double value)
{
    clear();

    m_type = Double;
    m_value._double = value;

    return true;
}

template<>
inline bool Variant::setValue(std::string string)
{
    clear();

    m_type = String;
    m_value.string = new std::string(string);

    return true;
}

template<>
inline bool Variant::setValue(const char *string)
{
    return setValue(std::string(string));
}

template<>
inline bool Variant::setValue(void *pointer)
{
    clear();

    m_type = Pointer;
    m_value.pointer = pointer;

    return true;
}

/// Returns the value of the variant in the type given by \c T.
template<typename T>
inline T Variant::value() const
{
    return 0;
}

template<>
inline bool Variant::value() const
{
    if(m_type == Bool){
        return m_value._bool;
    }
    else if(m_type == Int){
        return m_value._int != 0;
    }

    return false;
}

template<>
inline char Variant::value() const
{
    if(m_type == Int){
        return static_cast<char>(m_value._int);
    }
    else if(m_type == String && !m_value.string->empty()){
        return m_value.string->at(0);
    }

    return '\0';
}

template<>
inline short Variant::value() const
{
    if(m_type == Int){
        return static_cast<short>(m_value._int);
    }
    else if(m_type == String){
        return boost::lexical_cast<short>(*m_value.string);
    }

    return 0;
}

template<>
inline int Variant::value() const
{
    if(m_type == Int){
        return m_value._int;
    }
    else if(m_type == Long){
        return static_cast<int>(m_value._long);
    }
    else if(m_type == Bool){
        return static_cast<int>(m_value._bool);
    }
    else if(m_type == Float){
        return static_cast<int>(m_value._float);
    }
    else if(m_type == Double){
        return static_cast<int>(m_value._double);
    }
    else if(m_type == String){
        return boost::lexical_cast<int>(*m_value.string);
    }

    return 0;
}

template<>
inline long Variant::value() const
{
    if(m_type == Long){
        return m_value._long;
    }
    else if(m_type == Int){
        return static_cast<long>(m_value._int);
    }
    else if(m_type == String){
        return boost::lexical_cast<long>(*m_value.string);
    }

    return 0;
}

template<>
inline size_t Variant::value() const
{
    if(m_type == Long){
        return static_cast<size_t>(m_value._long);
    }
    else if(m_type == Int){
        return static_cast<size_t>(m_value._int);
    }
    else if(m_type == String){
        return boost::lexical_cast<size_t>(*m_value.string);
    }

    return 0;
}

template<>
inline float Variant::value() const
{
    if(m_type == Float){
        return m_value._float;
    }
    else if(m_type == Double){
        return static_cast<float>(m_value._double);
    }
    else if(m_type == Int){
        return static_cast<float>(m_value._int);
    }
    else if(m_type == Long){
        return static_cast<float>(m_value._long);
    }
    else if(m_type == String){
        return boost::lexical_cast<float>(*m_value.string);
    }

    return 0;
}

template<>
inline double Variant::value() const
{
    if(m_type == Double){
        return m_value._double;
    }
    else if(m_type == Float){
        return static_cast<double>(m_value._float);
    }
    else if(m_type == Int){
        return static_cast<double>(m_value._int);
    }
    else if(m_type == Long){
        return static_cast<double>(m_value._long);
    }
    else if(m_type == String){
        return boost::lexical_cast<double>(*m_value.string);
    }

    return 0;
}

template<>
inline void* Variant::value() const
{
    if(m_type == Pointer){
        return m_value.pointer;
    }

    return 0;
}

template<>
inline std::string Variant::value() const
{
    if(m_type == String){
        return *m_value.string;
    }
    else if(m_type == Int){
        return boost::lexical_cast<std::string>(m_value._int);
    }
    else if(m_type == Float){
        return boost::lexical_cast<std::string>(m_value._float);
    }
    else if(m_type == Double){
        return boost::lexical_cast<std::string>(m_value._double);
    }

    return std::string();
}

/// Clears the variant's data and sets the variant to null.
inline void Variant::clear()
{
    if(m_type == String){
        delete m_value.string;
        m_value.string = 0;
    }

    m_type = Null;
}

// --- Conversions --------------------------------------------------------- //
/// Returns the value of the variant as a \c bool.
inline bool Variant::toBool() const
{
    return value<bool>();
}

/// Returns the value of the variant as a \c char.
inline char Variant::toChar() const
{
    return value<char>();
}

/// Returns the value of the variant as an \c unsigned \c char.
inline unsigned char Variant::toUChar() const
{
    return value<unsigned char>();
}

/// Returns the value of the variant as a \c short.
inline short Variant::toShort() const
{
    return value<short>();
}

/// Returns the value of the variant as an \c unsigned \c short.
inline unsigned short Variant::toUShort() const
{
    return value<unsigned short>();
}

/// Returns the value of the variant as an \c int.
inline int Variant::toInt() const
{
    return value<int>();
}

/// Returns the value of the variant as an \c unsigned \c int.
inline unsigned int Variant::toUInt() const
{
    return value<unsigned int>();
}

/// Returns the value of the variant as a \c long.
inline long Variant::toLong() const
{
    return value<long>();
}

/// Returns the value of the variant as an \c unsigned \c long.
inline unsigned long Variant::toULong() const
{
    return value<unsigned long>();
}

/// Returns the value of the variant as a \p size_t.
inline size_t Variant::toSizeT() const
{
    return value<size_t>();
}

/// Returns the value of the variant as a \c float.
inline float Variant::toFloat() const
{
    return value<float>();
}

/// Returns the value of the variant as a \c double.
inline double Variant::toDouble() const
{
    return value<double>();
}

/// Returns the value of the variant as a \c Real.
inline Real Variant::toReal() const
{
    return value<Real>();
}

/// Returns the value of the variant as a pointer.
inline void* Variant::toPointer() const
{
    return value<void *>();
}

/// Returns the value of the variant as a string.
inline std::string Variant::toString() const
{
    return value<std::string>();
}

// --- Operators ----------------------------------------------------------- //
inline Variant& Variant::operator=(const Variant &variant)
{
    if(this != &variant){
        // clear previous data
        clear();

        // set new type
        m_type = variant.m_type;

        // set new value
        if(m_type == String){
            m_value.string = new std::string(variant.toString());
        }
        else if(m_type != Null){
            m_value = variant.m_value;
        }
    }

    return *this;
}

} // end chemkit namespace

#endif // CHEMKIT_VARIANT_INLINE_H
