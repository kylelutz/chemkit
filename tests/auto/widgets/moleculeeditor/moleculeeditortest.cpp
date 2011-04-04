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

#include "moleculeeditortest.h"

#include <chemkit/moleculeeditor.h>

void MoleculeEditorTest::basic()
{
    chemkit::MoleculeEditor editor;
    QVERIFY(editor.molecule() == 0);

    chemkit::Molecule molecule;
    chemkit::MoleculeEditor editor2(&molecule);
    QVERIFY(editor2.molecule() == &molecule);
}

void MoleculeEditorTest::setMolecule()
{
    chemkit::MoleculeEditor editor;
    QVERIFY(editor.molecule() == 0);

    chemkit::Molecule molecule;
    editor.setMolecule(&molecule);
    QVERIFY(editor.molecule() == &molecule);

    chemkit::Molecule molecule2;
    editor.setMolecule(&molecule2);
    QVERIFY(editor.molecule() == &molecule2);

    editor.setMolecule(0);
    QVERIFY(editor.molecule() == 0);
}

void MoleculeEditorTest::addAtom()
{
    chemkit::Molecule molecule;
    chemkit::MoleculeEditor editor(&molecule);
    QVERIFY(editor.molecule() == &molecule);

    chemkit::Atom *C1 = editor.addAtom(6);
    chemkit::Atom *C2 = editor.addAtom(6);
    chemkit::Atom *C3 = editor.addAtom(6);
    QCOMPARE(molecule.formula(), std::string("C3"));

    editor.undo();
    QCOMPARE(molecule.formula(), std::string("C2"));

    editor.redo();
    QCOMPARE(molecule.formula(), std::string("C3"));

    editor.removeAtom(C1);
    QCOMPARE(molecule.formula(), std::string("C2"));

    editor.undo();
    QCOMPARE(molecule.formula(), std::string("C3"));

    editor.redo();
    QCOMPARE(molecule.formula(), std::string("C2"));

    editor.removeAtom(C2);
    editor.removeAtom(C3);
    QCOMPARE(molecule.formula(), std::string());

    editor.undo();
    QCOMPARE(molecule.formula(), std::string("C"));

    editor.undo();
    QCOMPARE(molecule.formula(), std::string("C2"));
}

void MoleculeEditorTest::removeAtom()
{
}

void MoleculeEditorTest::setAtomAtomicNumber()
{
    chemkit::Molecule molecule;
    chemkit::Atom *atom = molecule.addAtom(6);
    chemkit::MoleculeEditor editor(&molecule);
    editor.setAtomAtomicNumber(atom, 1);
    QCOMPARE(atom->atomicNumber(), 1);

    editor.undo();
    QCOMPARE(atom->atomicNumber(), 6);

    editor.redo();
    QCOMPARE(atom->atomicNumber(), 1);

    editor.setAtomAtomicNumber(atom, 2);
    QCOMPARE(atom->atomicNumber(), 2);

    editor.setAtomAtomicNumber(atom, 3);
    QCOMPARE(atom->atomicNumber(), 3);

    editor.undo();
    QCOMPARE(atom->atomicNumber(), 2);

    editor.undo();
    QCOMPARE(atom->atomicNumber(), 1);
}

void MoleculeEditorTest::setAtomPosition()
{
}

void MoleculeEditorTest::addBond()
{
    chemkit::Molecule molecule;
    chemkit::MoleculeEditor editor(&molecule);
    QVERIFY(editor.molecule() == &molecule);

    chemkit::Atom *C1 = editor.addAtom(6);
    chemkit::Atom *C2 = editor.addAtom(6);
    chemkit::Atom *C3 = editor.addAtom(6);
    QCOMPARE(molecule.formula(), std::string("C3"));

    chemkit::Bond *C1_C2 = editor.addBond(C1, C2);
    QCOMPARE(molecule.bondCount(), 1);

    chemkit::Bond *C2_C3 = editor.addBond(C2, C3, chemkit::Bond::Double);
    QCOMPARE(C2_C3->order(), 2);
    QCOMPARE(molecule.bondCount(), 2);

    editor.removeBond(C1_C2);
    QCOMPARE(molecule.bondCount(), 1);

    editor.undo();
    QCOMPARE(molecule.bondCount(), 2);

    editor.redo();
    QCOMPARE(molecule.bondCount(), 1);

    editor.removeBond(C2_C3);
    QCOMPARE(molecule.bondCount(), 0);

    editor.undo();
    QCOMPARE(molecule.bondCount(), 1);
    QCOMPARE(molecule.bonds()[0]->order(), 2);
}

void MoleculeEditorTest::removeBond()
{
    chemkit::Molecule molecule;
    chemkit::MoleculeEditor editor(&molecule);
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(C1, C2);
    QCOMPARE(molecule.bondCount(), 1);

    editor.removeBond(bond);
    QCOMPARE(molecule.bondCount(), 0);
    QCOMPARE(C1->isBondedTo(C2), false);

    editor.undo();
    QCOMPARE(molecule.bondCount(), 1);
    QCOMPARE(C1->isBondedTo(C2), true);

    editor.redo();
    QCOMPARE(molecule.bondCount(), 0);
    QCOMPARE(C1->isBondedTo(C2), false);
}

void MoleculeEditorTest::setBondOrder()
{
    chemkit::Molecule molecule;
    chemkit::MoleculeEditor editor(&molecule);
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Bond *bond = molecule.addBond(C1, C2);
    QCOMPARE(bond->order(), 1);

    editor.setBondOrder(bond, 2);
    QCOMPARE(bond->order(), 2);

    editor.undo();
    QCOMPARE(bond->order(), 1);

    editor.redo();
    QCOMPARE(bond->order(), 2);

    editor.setBondOrder(bond, 3);
    editor.setBondOrder(bond, 2);
    QCOMPARE(bond->order(), 2);

    editor.undo();
    QCOMPARE(bond->order(), 3);

    editor.undo();
    QCOMPARE(bond->order(), 2);
}

void MoleculeEditorTest::copy()
{
    chemkit::Molecule molecule;
    chemkit::Atom *C1 = molecule.addAtom("C");
    chemkit::Atom *C2 = molecule.addAtom("C");
    chemkit::Atom *O3 = molecule.addAtom("O");
    molecule.addBond(C1, C2);
    molecule.addBond(C2, O3, 2);
    QCOMPARE(molecule.formula(), std::string("C2O"));
    QCOMPARE(molecule.bondCount(), 2);

    chemkit::MoleculeEditor editor(&molecule);
    editor.copy(molecule.atoms());
    QCOMPARE(editor.copyBuffer().size(), 3);

    editor.paste();
    QCOMPARE(molecule.formula(), std::string("C4O2"));
    QCOMPARE(molecule.bondCount(), 4);
    foreach(chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            QCOMPARE(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Double), true);
        }
        else if(atom->is(chemkit::Atom::Carbon)){
            QCOMPARE(atom->isBondedTo(chemkit::Atom::Carbon, chemkit::Bond::Single), true);
        }
    }

    editor.undo();
    QCOMPARE(molecule.formula(), std::string("C2O"));
    QCOMPARE(molecule.bondCount(), 2);

    editor.redo();
    QCOMPARE(molecule.formula(), std::string("C4O2"));
    QCOMPARE(molecule.bondCount(), 4);

    QList<chemkit::Atom *> oxygens;
    foreach(chemkit::Atom *atom, molecule.atoms()){
        if(atom->is(chemkit::Atom::Oxygen)){
            oxygens.append(atom);
        }
    }
    editor.copy(oxygens);
    QCOMPARE(editor.copyBuffer().size(), 2);

    QList<chemkit::Atom *> newAtoms = editor.paste();
    QCOMPARE(newAtoms.size(), 2);
    foreach(chemkit::Atom *atom, newAtoms){
        QCOMPARE(atom->is(chemkit::Atom::Oxygen), true);
    }
    QCOMPARE(molecule.formula(), std::string("C4O4"));
    QCOMPARE(molecule.bondCount(), 4);
}

void MoleculeEditorTest::clearCopyBuffer()
{
    chemkit::MoleculeEditor editor;
    QCOMPARE(editor.copyBuffer().size(), 0);
    editor.clearCopyBuffer();
    QCOMPARE(editor.copyBuffer().size(), 0);

    chemkit::Molecule molecule;
    molecule.addAtom("H");
    molecule.addAtom("H");
    editor.setMolecule(&molecule);
    editor.copy(molecule.atoms());
    QCOMPARE(editor.copyBuffer().size(), 2);
    editor.clearCopyBuffer();
    QCOMPARE(editor.copyBuffer().size(), 0);
}

QTEST_APPLESS_MAIN(MoleculeEditorTest)
