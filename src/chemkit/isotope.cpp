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

#include "isotope.h"

namespace chemkit {

// === Isotope ============================================================= //
/// \class Isotope isotope.h chemkit/isotope.h
/// \ingroup chemkit
/// \brief The Isotope class represents an isotope.

// --- Construction and Destruction ---------------------------------------- //
/// Create a new, invalid isotope.
Isotope::Isotope()
{
    m_neutronCount = 0;
}

/// Create a new isotope with \p element and an identical number of
/// protons and neutrons.
Isotope::Isotope(const Element &element)
    : m_element(element)
{
    m_neutronCount = element.atomicNumber();
}

/// Create a new isotope with \p element and \p massNumber.
Isotope::Isotope(const Element &element, size_t massNumber)
    : m_element(element)
{
    setMassNumber(massNumber);
}

// --- Properties ---------------------------------------------------------- //
/// Sets the element for the isotope to \p element.
void Isotope::setElement(const Element &element)
{
    m_element = element;
}

/// Returns the element for the isotope.
Element Isotope::element() const
{
    return m_element;
}

/// Sets the number of protons in the isotope to \p count.
void Isotope::setProtonCount(size_t count)
{
    m_element.setAtomicNumber(count);
}

/// Returns the number of protons in the isotope.
size_t Isotope::protonCount() const
{
    return m_element.atomicNumber();
}

/// Sets the number of neutrons to \p count.
void Isotope::setNeutronCount(size_t count)
{
    m_neutronCount = count;
}

/// Returns the number of neutrons in the isotope.
size_t Isotope::neutronCount() const
{
    return m_neutronCount;
}

/// Sets the atomic number for the isotop to \p number.
void Isotope::setAtomicNumber(AtomicNumberType number)
{
    m_element.setAtomicNumber(number);
}

/// Returns the atomic number for the isotope.
Isotope::AtomicNumberType Isotope::atomicNumber() const
{
    return m_element.atomicNumber();
}

/// Sets the mass number of the isotope to \p number.
void Isotope::setMassNumber(MassNumberType number)
{
    m_neutronCount = static_cast<unsigned char>(number - protonCount());
}

/// Return the mass number for the isotope. This is equal to the
/// number of protons plus the number of neutrons.
Isotope::MassNumberType Isotope::massNumber() const
{
    return static_cast<MassNumberType>(protonCount() + neutronCount());
}

} // end chemkit namespace
