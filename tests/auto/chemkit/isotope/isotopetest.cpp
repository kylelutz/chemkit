/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "isotopetest.h"

#include <chemkit/isotope.h>

void IsotopeTest::constructor()
{
    chemkit::Isotope isotope;
    QCOMPARE(isotope.element(), chemkit::Element());

    chemkit::Isotope C12("C");
    QCOMPARE(C12.element(), chemkit::Element("C"));
    QCOMPARE(C12.protonCount(), size_t(6));
    QCOMPARE(C12.neutronCount(), size_t(6));

    chemkit::Isotope U238("U", 238);
    QCOMPARE(U238.atomicNumber(), chemkit::Isotope::AtomicNumberType(92));
    QCOMPARE(U238.massNumber(), chemkit::Isotope::MassNumberType(238));

    chemkit::Isotope C12_copy(C12);
    QCOMPARE(C12_copy.element(), chemkit::Element("C"));
    QCOMPARE(C12_copy.protonCount(), size_t(6));
    QCOMPARE(C12_copy.neutronCount(), size_t(6));
}

void IsotopeTest::element()
{
    chemkit::Isotope isotope;
    QCOMPARE(isotope.element(), chemkit::Element());

    isotope.setElement("N");
    QCOMPARE(isotope.element(), chemkit::Element("N"));

    isotope.setElement(chemkit::Element());
    QCOMPARE(isotope.element(), chemkit::Element());
}

void IsotopeTest::atomicNumber()
{
    chemkit::Isotope isotope;
    QCOMPARE(isotope.atomicNumber(), chemkit::Isotope::AtomicNumberType(0));

    isotope.setAtomicNumber(16);
    QCOMPARE(isotope.atomicNumber(), chemkit::Isotope::AtomicNumberType(16));
}

void IsotopeTest::massNumber()
{
    chemkit::Isotope isotope;
    QCOMPARE(isotope.massNumber(), chemkit::Isotope::MassNumberType(0));

    isotope.setMassNumber(10);
    QCOMPARE(isotope.massNumber(), chemkit::Isotope::MassNumberType(10));

    isotope.setElement("B");
    QCOMPARE(isotope.atomicNumber(), chemkit::Isotope::AtomicNumberType(5));
    QCOMPARE(isotope.massNumber(), chemkit::Isotope::MassNumberType(15));
}

void IsotopeTest::neutronCount()
{
    chemkit::Isotope isotope;
    QCOMPARE(isotope.neutronCount(), size_t(0));

    isotope.setNeutronCount(20);
    QCOMPARE(isotope.neutronCount(), size_t(20));
}

QTEST_APPLESS_MAIN(IsotopeTest)
