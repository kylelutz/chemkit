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

#ifndef MOLECULARSURFACETEST_H
#define MOLECULARSURFACETEST_H

#include <QtTest>

class MolecularSurfaceTest : public QObject
{
    Q_OBJECT

    private slots:
        // function tests
        void molecule();
        void probeRadius();
        void surfaceType();

        // molecule tests
        void hydrogen();
        void water();
        void serine();
        void guanine();
        void methane();
        void ethanol();
        void adenosine();
        void lysozyme();
        void cytochrome();
        void toxin();
        void hydrolase();
        void hemoglobin();
        void dna();
        void ribozyme();
        void ubiqutin();
};

#endif // MOLECULARSURFACETEST_H
