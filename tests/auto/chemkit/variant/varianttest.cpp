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

#include "varianttest.h"

#include <chemkit/variant.h>

void VariantTest::isNull()
{
    chemkit::Variant variant;
    QCOMPARE(variant.isNull(), true);

    variant.setValue(7);
    QCOMPARE(variant.isNull(), false);
}

void VariantTest::clear()
{
    chemkit::Variant variant(62);
    QCOMPARE(variant.isNull(), false);

    variant.clear();
    QCOMPARE(variant.isNull(), true);

    variant.setValue('f');
    QCOMPARE(variant.isNull(), false);

    variant.clear();
    QCOMPARE(variant.isNull(), true);
}

void VariantTest::toBool()
{
    chemkit::Variant variant(false);
    QCOMPARE(variant.toBool(), false);

    variant.setValue(true);
    QCOMPARE(variant.toBool(), true);

    variant.setValue(0);
    QCOMPARE(variant.toBool(), false);

    variant.setValue(1);
    QCOMPARE(variant.toBool(), true);

    variant.setValue(-5);
    QCOMPARE(variant.toBool(), true);
}

void VariantTest::toChar()
{
    chemkit::Variant variant('c');
    QCOMPARE(variant.toChar(), 'c');

    variant.setValue("hello");
    QCOMPARE(variant.toChar(), 'h');
}

void VariantTest::toShort()
{
    chemkit::Variant variant(short(4));
    QCOMPARE(variant.toShort(), short(4));
}

void VariantTest::toInt()
{
    chemkit::Variant variant(12);
    QCOMPARE(variant.toInt(), int(12));

    variant.setValue(-23);
    QCOMPARE(variant.toInt(), int(-23));

    variant.setValue("42");
    QCOMPARE(variant.toInt(), int(42));

    variant.setValue(true);
    QCOMPARE(variant.toInt(), int(1));

    variant.setValue(false);
    QCOMPARE(variant.toInt(), int(0));
}

void VariantTest::toLong()
{
    chemkit::Variant variant(192L);
    QCOMPARE(variant.toLong(), 192L);

    variant.setValue(7);
    QCOMPARE(variant.toLong(), 7L);

    variant.setValue("562");
    QCOMPARE(variant.toLong(), 562L);
}

void VariantTest::toSizeT()
{
    chemkit::Variant variant(size_t(12345));
    QCOMPARE(variant.toSizeT(), size_t(12345));

    variant.setValue(size_t(-54321));
    QCOMPARE(variant.toSizeT(), size_t(-54321));

    variant.setValue("98765");
    QCOMPARE(variant.toSizeT(), size_t(98765));
}

void VariantTest::toFloat()
{
    chemkit::Variant variant(12.3f);
    QCOMPARE(variant.toFloat(), 12.3f);
}

void VariantTest::toDouble()
{
    chemkit::Variant variant(3.14);
    QCOMPARE(variant.toDouble(), 3.14);
}

void VariantTest::toPointer()
{
    int value;
    void *pointer = &value;
    chemkit::Variant variant(pointer);
    QCOMPARE(variant.toPointer(), pointer);
}

void VariantTest::toString()
{
    chemkit::Variant variant("hello");
    QCOMPARE(variant.toString(), std::string("hello"));

    variant.setValue(12);
    QCOMPARE(variant.toString(), std::string("12"));

    variant.setValue(std::string("hello2"));
    QCOMPARE(variant.toString(), std::string("hello2"));
}

QTEST_APPLESS_MAIN(VariantTest)
