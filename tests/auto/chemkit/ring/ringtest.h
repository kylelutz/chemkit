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

#ifndef RINGTEST_H
#define RINGTEST_H

#include <QtTest>

#include <chemkit/ring.h>
#include <chemkit/molecule.h>

class RingTest : public QObject
{
    Q_OBJECT

    private:
        chemkit::Molecule *benzene;
        chemkit::Molecule *furan;
        chemkit::Molecule *cyclohexane;
        chemkit::Molecule *cyclopropane;

        chemkit::Ring *benzeneRing;
        chemkit::Ring *furanRing;
        chemkit::Ring *cyclohexaneRing;
        chemkit::Ring *cyclopropaneRing;

    private slots:
        void initTestCase();
        void molecule();
        void atoms();
        void atomCount();
        void size();
        void bonds();
        void bondCount();
        void root();
        void position();
        void contains();
        void heteroatomCount();
        void isHeterocycle();
        void isAromatic();
        void cleanupTestCase();
};

#endif // RINGTEST_H
