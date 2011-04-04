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

#ifndef SUBSTRUCTURESEARCHTEST_H
#define SUBSTRUCTURESEARCHTEST_H

#include <QtTest>

#include <chemkit/molecule.h>

class SubstructureSearchTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();
        void cleanupTestCase();

        void benzene();
        void butane();
        void cyclopropane();
        void ethane();
        void ethanol();
        void indole();
        void methane();
        void methanol();
        void phenol();
        void propane();

        void protein();

    private:
        chemkit::Molecule *m_benzene;
        chemkit::Molecule *m_butane;
        chemkit::Molecule *m_cyclopropane;
        chemkit::Molecule *m_ethane;
        chemkit::Molecule *m_ethanol;
        chemkit::Molecule *m_indole;
        chemkit::Molecule *m_methane;
        chemkit::Molecule *m_methanol;
        chemkit::Molecule *m_phenol;
        chemkit::Molecule *m_propane;
};

#endif // SUBSTRUCTURESEARCHTEST_H
