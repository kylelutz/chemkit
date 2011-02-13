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

#include <QtTest>

#include <chemkit/chemkit.h>
#include <chemkit/residue.h>
#include <chemkit/molecule.h>

class ResidueTest : public QObject
{
    Q_OBJECT

    private slots:
        void molecule();
        void atomType();
};

void ResidueTest::molecule()
{
    chemkit::Molecule molecule;
    chemkit::Residue *residue = new chemkit::Residue(&molecule);
    QVERIFY(residue->molecule() == &molecule);
}

void ResidueTest::atomType()
{
    chemkit::Molecule molecule;
    chemkit::Residue *residue = new chemkit::Residue(&molecule);

    chemkit::Atom *c1 = molecule.addAtom("C");
    chemkit::Atom *c2 = molecule.addAtom("C");
    residue->addAtom(c1);
    residue->addAtom(c2);
    QCOMPARE(residue->atomType(c1), QString());
    QCOMPARE(residue->atomType(c2), QString());

    residue->setAtomType(c1, "C1");
    QCOMPARE(residue->atomType(c1), QString("C1"));
    QVERIFY(residue->atom("C1") == c1);
    QVERIFY(residue->atom("C2") == 0);

    residue->setAtomType(c2, "C2");
    QCOMPARE(residue->atomType(c2), QString("C2"));
    QVERIFY(residue->atom("C2") == c2);
}

QTEST_APPLESS_MAIN(ResidueTest)
#include "residuetest.moc"
