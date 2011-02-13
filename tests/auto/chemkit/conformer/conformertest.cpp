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

#include <chemkit/molecule.h>
#include <chemkit/conformer.h>

class ConformerTest : public QObject
{
    Q_OBJECT

    private slots:
        void basic();
        void atomPosition();
};

void ConformerTest::basic()
{
    chemkit::Molecule molecule;
    chemkit::Conformer *conformer = molecule.conformer();
    QVERIFY(conformer->molecule() == &molecule);

    conformer = molecule.addConformer();
    QVERIFY(conformer->molecule() == &molecule);
}

void ConformerTest::atomPosition()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom("C");
    atom->setPosition(1, 2, 3);
    QCOMPARE(atom->position(), chemkit::Point(1, 2, 3));

    chemkit::Conformer *conformer = molecule.addConformer();
    QVERIFY(conformer != 0);

    conformer->setPosition(atom, chemkit::Point(3, 2, 1));
    QCOMPARE(conformer->position(atom), chemkit::Point(3, 2, 1));
    QCOMPARE(atom->position(), chemkit::Point(1, 2, 3));

    chemkit::Conformer *originalConformer = molecule.conformer();
    QVERIFY(originalConformer != 0);

    molecule.setConformer(conformer);
    QCOMPARE(atom->position(), chemkit::Point(3, 2, 1));
    QCOMPARE(originalConformer->position(atom), chemkit::Point(1, 2, 3));
}

QTEST_APPLESS_MAIN(ConformerTest)
#include "conformertest.moc"
