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

#ifndef RINGPERCEPTIONTEST_H
#define RINGPERCEPTIONTEST_H

#include <QtTest>

#include <chemkit/molecule.h>

class RingPerceptionTest : public QObject
{
    Q_OBJECT

    private:
        static void addHydrogens(chemkit::Molecule *molecule);

    private slots:
        void anthracene();
        void anthraquinone();
        void arsole();
        void benzene();
        void benzimidazole();
        void benzobicyclooctane();
        void benzonorborene();
        void bicyclooctane();
        void biotin();
        void biphenylene();
        void cubane();
        void cyclobutane();
        void cyclodecane();
        void cycloheptane();
        void cyclohexane();
        void cyclononane();
        void cyclooctane();
        void cyclooctatetraene();
        void cyclopentane();
        void cyclopropane();
        void decalin();
        void furan();
        void indole();
        void imidazole();
        void ladderane();
        void naphthalene();
        void norbornane();
        void oxazole();
        void oxirane();
        void porphin();
        void pyrazole();
        void pyrene();
        void pyridine();
        void pyrrole();
        void quinoxaline();
        void tetralin();
        void thiophene();
        void tricyclohexane();
        void tricyclooctane();
        void uracil();
        void vigtua();
};

#endif // RINGPERCEPTIONTEST_H
