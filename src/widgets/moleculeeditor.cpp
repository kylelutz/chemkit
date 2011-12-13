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

#include "moleculeeditor.h"

#include <cassert>

#include <QHash>
#include <QVector>
#include <QUndoCommand>

#include <chemkit/atom.h>
#include <chemkit/bond.h>

namespace chemkit {

namespace {

// === MoleculeEditorCommand =============================================== //
class MoleculeEditorCommand : public QUndoCommand
{
public:
    MoleculeEditorCommand(MoleculeEditor *editor);
    virtual ~MoleculeEditorCommand();

    MoleculeEditor* editor() const { return m_editor; }
    Molecule* molecule() const { return m_editor->molecule(); }

private:
    MoleculeEditor *m_editor;
};

MoleculeEditorCommand::MoleculeEditorCommand(MoleculeEditor *editor)
    : QUndoCommand()
{
    m_editor = editor;
}

MoleculeEditorCommand::~MoleculeEditorCommand()
{
}

// === AddAtomCommand ====================================================== //
class AddAtomCommand : public MoleculeEditorCommand
{
public:
    AddAtomCommand(MoleculeEditor *editor, const Element &element);

    void undo();
    void redo();

    Atom* atom() const { return m_atom; }

private:
    Element m_element;
    Atom *m_atom;
    int m_atomId;
};

AddAtomCommand::AddAtomCommand(MoleculeEditor *editor, const Element &element)
    : MoleculeEditorCommand(editor),
      m_element(element)
{
    m_atomId = 0;
}

void AddAtomCommand::undo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    molecule()->removeAtom(atom);
}

void AddAtomCommand::redo()
{
    m_atom = molecule()->addAtom(m_element);

    if(m_atomId){
        editor()->setId(m_atom, m_atomId);
    }
    else{
        m_atomId = editor()->id(m_atom);
    }
}

// === RemoveAtomCommand =================================================== //
class RemoveAtomCommand : public MoleculeEditorCommand
{
public:
    RemoveAtomCommand(MoleculeEditor *editor, Atom *atom);

    void undo();
    void redo();

private:
    int m_atomId;
    int m_atomicNumber;
    Point3 m_position;
};

RemoveAtomCommand::RemoveAtomCommand(MoleculeEditor *editor, Atom *atom)
    : MoleculeEditorCommand(editor)
{
    m_atomId = editor->id(atom);
    m_atomicNumber = atom->atomicNumber();
    m_position = atom->position();
}

void RemoveAtomCommand::undo()
{
    Atom *atom = molecule()->addAtom(m_atomicNumber);
    assert(atom);

    atom->setPosition(m_position);
    editor()->setId(atom, m_atomId);
}

void RemoveAtomCommand::redo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    molecule()->removeAtom(atom);
}

// === SetAtomAtomicNumberCommand ========================================== //
class SetAtomAtomicNumberCommand : public MoleculeEditorCommand
{
public:
    SetAtomAtomicNumberCommand(MoleculeEditor *editor, Atom *atom, int atomicNumber);

    void undo();
    void redo();

private:
    int m_atomId;
    int m_initialAtomicNumber;
    int m_finalAtomicNumber;
};

SetAtomAtomicNumberCommand::SetAtomAtomicNumberCommand(MoleculeEditor *editor, Atom *atom, int atomicNumber)
    : MoleculeEditorCommand(editor)
{
    m_atomId = editor->id(atom);
    m_initialAtomicNumber = atom->atomicNumber();
    m_finalAtomicNumber = atomicNumber;
}

void SetAtomAtomicNumberCommand::undo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    atom->setAtomicNumber(m_initialAtomicNumber);
}

void SetAtomAtomicNumberCommand::redo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    atom->setAtomicNumber(m_finalAtomicNumber);
}

// === SetAtomPositionCommand ============================================== //
class SetAtomPositionCommand : public MoleculeEditorCommand
{
public:
    SetAtomPositionCommand(MoleculeEditor *editor, Atom *atom, const Point3 &position);

    void undo();
    void redo();

private:
    int m_atomId;
    Point3 m_initialPosition;
    Point3 m_finalPosition;
};

SetAtomPositionCommand::SetAtomPositionCommand(MoleculeEditor *editor, Atom *atom, const Point3 &position)
    : MoleculeEditorCommand(editor)
{
    m_atomId = editor->id(atom);
    m_initialPosition = atom->position();
    m_finalPosition = position;
}

void SetAtomPositionCommand::undo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    atom->setPosition(m_initialPosition);
}

void SetAtomPositionCommand::redo()
{
    Atom *atom = editor()->atom(m_atomId);
    assert(atom);

    atom->setPosition(m_finalPosition);
}

// === AddBondCommand ====================================================== //
class AddBondCommand : public MoleculeEditorCommand
{
public:
    AddBondCommand(MoleculeEditor *editor, Atom *a, Atom *b);

    void undo();
    void redo();

    Bond* bond() const { return m_bond; }

private:
    int m_atomId1;
    int m_atomId2;
    Bond *m_bond;
};

AddBondCommand::AddBondCommand(MoleculeEditor *editor, Atom *a, Atom *b)
    : MoleculeEditorCommand(editor)
{
    m_atomId1 = editor->id(a);
    m_atomId2 = editor->id(b);
}

void AddBondCommand::undo()
{
    Bond *bond = editor()->bond(m_atomId1, m_atomId2);
    assert(bond);

    molecule()->removeBond(bond);
}

void AddBondCommand::redo()
{
    Atom *atom1 = editor()->atom(m_atomId1);
    Atom *atom2 = editor()->atom(m_atomId2);
    assert(atom1);
    assert(atom2);

    m_bond = molecule()->addBond(atom1, atom2);
}

// === RemoveBondCommand =================================================== //
class RemoveBondCommand : public MoleculeEditorCommand
{
public:
    RemoveBondCommand(MoleculeEditor *editor, Bond *bond);

    void undo();
    void redo();

private:
    int m_atomId1;
    int m_atomId2;
    int m_bondOrder;
};

RemoveBondCommand::RemoveBondCommand(MoleculeEditor *editor, Bond *bond)
    : MoleculeEditorCommand(editor)
{
    m_atomId1 = editor->id(bond->atom1());
    m_atomId2 = editor->id(bond->atom2());
    m_bondOrder = bond->order();
}

void RemoveBondCommand::undo()
{
    Atom *atom1 = editor()->atom(m_atomId1);
    Atom *atom2 = editor()->atom(m_atomId2);
    assert(!atom1->isBondedTo(atom2));

    molecule()->addBond(atom1, atom2, m_bondOrder);
}

void RemoveBondCommand::redo()
{
    Bond *bond = editor()->bond(m_atomId1, m_atomId2);
    assert(bond != 0);

    molecule()->removeBond(bond);
}

// === SetBondOrderCommand ================================================= //
class SetBondOrderCommand : public MoleculeEditorCommand
{
public:
    SetBondOrderCommand(MoleculeEditor *editor, Bond *bond, int order);

    void undo();
    void redo();

private:
    int m_atomId1;
    int m_atomId2;
    int m_initialOrder;
    int m_finalOrder;
};

SetBondOrderCommand::SetBondOrderCommand(MoleculeEditor *editor, Bond *bond, int order)
    : MoleculeEditorCommand(editor)
{
    m_atomId1 = editor->id(bond->atom1());
    m_atomId2 = editor->id(bond->atom2());
    m_initialOrder = bond->order();
    m_finalOrder = order;
}

void SetBondOrderCommand::undo()
{
    Bond *bond = editor()->bond(m_atomId1, m_atomId2);
    assert(bond != 0);

    bond->setOrder(m_initialOrder);
}

void SetBondOrderCommand::redo()
{
    Bond *bond = editor()->bond(m_atomId1, m_atomId2);
    assert(bond != 0);

    bond->setOrder(m_finalOrder);
}

} // end anonymous namespace

// === MoleculeEditorPrivate =============================================== //
class MoleculeEditorPrivate
{
public:
    Molecule *molecule;
    bool inEdit;
    QUndoStack undoStack;
    QHash<int, Atom *> atomIds;
    QList<Atom *> copyBuffer;
    Molecule *cutMolecule;
};

// === MoleculeEditor ====================================================== //
/// \class MoleculeEditor moleculeeditor.h chemkit/moleculeeditor.h
/// \ingroup chemkit-widgets
/// \brief The MoleculeEditor class provides editing functions for
///        molecules.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new molecule editor object for \p molecule.
MoleculeEditor::MoleculeEditor(Molecule *molecule)
    : QObject(),
      d(new MoleculeEditorPrivate)
{
    d->molecule = molecule;
    d->inEdit = false;
    d->cutMolecule = new Molecule;

    connect(&d->undoStack, SIGNAL(canUndoChanged(bool)), SLOT(canUndoChangedSlot(bool)));
    connect(&d->undoStack, SIGNAL(canRedoChanged(bool)), SLOT(canRedoChangedSlot(bool)));
}

/// Destroys the molecule editor object.
MoleculeEditor::~MoleculeEditor()
{
    delete d->cutMolecule;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the molecule to edit.
void MoleculeEditor::setMolecule(Molecule *molecule)
{
    if(molecule == d->molecule){
        return;
    }

    clearUndoStack();

    d->molecule = molecule;
}

/// Returns the molecule that is being edited.
Molecule* MoleculeEditor::molecule() const
{
    return d->molecule;
}

// --- Editing ------------------------------------------------------------- //
/// Reverts the last change the occured.
void MoleculeEditor::undo()
{
    d->undoStack.undo();
}

/// Returns \c true if it is possible to undo a change.
bool MoleculeEditor::canUndo() const
{
    return d->undoStack.canUndo();
}

/// Redoes the last action that was reverted by undo().
void MoleculeEditor::redo()
{
    d->undoStack.redo();
}

/// Returns \c true if it is possible to redo a change.
bool MoleculeEditor::canRedo() const
{
    return d->undoStack.canRedo();
}

/// Clears all undo actions.
void MoleculeEditor::clearUndoStack()
{
    d->undoStack.clear();
}

/// Starts an edit action. All modifications performed between
/// beginEdit() and endEdit() are grouped into a single action.
void MoleculeEditor::beginEdit()
{
    d->undoStack.beginMacro(QString());
    d->inEdit = true;
}

/// Ends an edit action.
void MoleculeEditor::endEdit()
{
    d->undoStack.endMacro();
    d->inEdit = false;
}

/// Returns \c true if the editor is in an edit action.
bool MoleculeEditor::isInEdit() const
{
    return d->inEdit;
}

/// Cuts each atom in \p atoms from the molecule.
void MoleculeEditor::cut(const QList<Atom *> &atoms)
{
    d->cutMolecule->clear();

    QHash<const Atom *, Atom *> cutAtomMap;

    foreach(Atom *atom, atoms){
        Atom *cutAtom = d->cutMolecule->addAtomCopy(atom);
        cutAtomMap[atom] = cutAtom;
    }

    for(int i = 0; i < atoms.size(); i++){
        for(int j = i + 1; j < atoms.size(); j++){
            Bond *bond = atoms[i]->bondTo(atoms[j]);
            if(bond)
                d->cutMolecule->addBond(cutAtomMap[atoms[i]], cutAtomMap[atoms[j]], bond->order());
        }
    }

    bool wasInEdit = isInEdit();

    if(!isInEdit()){
        beginEdit();
    }

    foreach(Atom *atom, atoms){
        removeAtom(atom);
    }

    if(!wasInEdit){
        endEdit();
    }

    d->copyBuffer = QVector<Atom *>::fromStdVector(d->cutMolecule->atoms()).toList();

    emit canPasteChanged(true);
}

/// Copies each atom in \p atoms.
void MoleculeEditor::copy(const QList<Atom *> &atoms)
{
    d->copyBuffer = atoms;

    emit canPasteChanged(true);
}

/// Paste the atoms from the copy buffer.
QList<Atom *> MoleculeEditor::paste()
{
    bool wasInEdit = isInEdit();

    if(!isInEdit()){
        beginEdit();
    }

    QHash<const Atom *, Atom *> oldToNew;

    foreach(const Atom *atom, d->copyBuffer){
        Atom *newAtom = addAtomCopy(atom);
        oldToNew[atom] = newAtom;
    }

    for(int i = 0; i < d->copyBuffer.size(); i++){
        Atom *oldAtom1 = d->copyBuffer[i];
        for(int j = i + 1; j < d->copyBuffer.size(); j++){
            Atom *oldAtom2 = d->copyBuffer[j];
            Bond *bond = oldAtom1->bondTo(oldAtom2);
            if(bond){
                addBond(oldToNew[oldAtom1], oldToNew[oldAtom2], bond->order());
            }
        }
    }

    if(!wasInEdit){
        endEdit();
    }

    return oldToNew.values();
}

/// Returns \c true if it is possible to paste atoms.
bool MoleculeEditor::canPaste() const
{
    return !copyBuffer().isEmpty();
}

/// Returns a list of atoms in the copy buffer.
QList<Atom *> MoleculeEditor::copyBuffer() const
{
    return d->copyBuffer;
}

/// Clears all atoms from the copy buffer.
void MoleculeEditor::clearCopyBuffer()
{
    d->copyBuffer.clear();

    emit canPasteChanged(false);
}

// --- Modification -------------------------------------------------------- //
/// Adds a new atom to the molecule.
///
/// \see Molecule::addAtom()
Atom* MoleculeEditor::addAtom(const Element &element)
{
    AddAtomCommand *command = new AddAtomCommand(this, element);
    d->undoStack.push(command);
    return command->atom();
}

/// Adds a copy of \p atom to the molecule.
///
/// \see Molecule::addAtomCopy()
Atom* MoleculeEditor::addAtomCopy(const Atom *atom)
{
    bool wasInEdit = isInEdit();

    if(!isInEdit()){
        beginEdit();
    }

    Atom *newAtom = addAtom(atom->atomicNumber());
    setAtomPosition(newAtom, atom->position());

    if(!wasInEdit){
        endEdit();
    }

    return newAtom;
}

/// Removes \p atom from the molecule.
///
/// \see Molecule::removeAtom()
void MoleculeEditor::removeAtom(Atom *atom)
{
    bool wasInEdit = isInEdit();

    if(!isInEdit()){
        beginEdit();
    }

    foreach(Bond *bond, atom->bonds()){
        removeBond(bond);
    }

    RemoveAtomCommand *command = new RemoveAtomCommand(this, atom);
    d->undoStack.push(command);

    if(!wasInEdit){
        endEdit();
    }
}

/// Sets the atomic number of \p atom to \p atomicNumber.
///
/// \see Atom::setAtomicNumber()
void MoleculeEditor::setAtomAtomicNumber(Atom *atom, int atomicNumber)
{
    SetAtomAtomicNumberCommand *command = new SetAtomAtomicNumberCommand(this, atom, atomicNumber);
    d->undoStack.push(command);
}

/// Sets the position of \p atom to \p position.
///
/// \see Atom::setPosition()
void MoleculeEditor::setAtomPosition(Atom *atom, const Point3 &position)
{
    SetAtomPositionCommand *command = new SetAtomPositionCommand(this, atom, position);
    d->undoStack.push(command);
}

/// Adds a bond between atoms \p a and \p b with \p order.
///
/// \see Molecule::addBond()
Bond* MoleculeEditor::addBond(Atom *a, Atom *b, int order)
{
    bool wasInEdit = isInEdit();

    if(!isInEdit()){
        beginEdit();
    }

    AddBondCommand *command = new AddBondCommand(this, a, b);
    d->undoStack.push(command);

    Bond *bond = command->bond();
    setBondOrder(bond, order);

    if(!wasInEdit){
        endEdit();
    }

    return bond;
}

/// Removes \p bond from the molecule.
///
/// \see Molecule::removeBond()
void MoleculeEditor::removeBond(Bond *bond)
{
    RemoveBondCommand *command = new RemoveBondCommand(this, bond);
    d->undoStack.push(command);
}

/// Sets the bond order for \p bond.
///
/// \see Bond::setOrder()
void MoleculeEditor::setBondOrder(Bond *bond, int order)
{
    SetBondOrderCommand *command = new SetBondOrderCommand(this, bond, order);
    d->undoStack.push(command);
}

// --- Internal Methods ---------------------------------------------------- //
Atom* MoleculeEditor::atom(int id)
{
    assert(d->atomIds.contains(id));

    return d->atomIds[id];
}

Bond* MoleculeEditor::bond(int id1, int id2)
{
    return atom(id1)->bondTo(atom(id2));
}

int MoleculeEditor::id(Atom *atom)
{
    int id = d->atomIds.key(atom, 0);

    if(!id){
        for(id = 1;;id++){
            if(!d->atomIds.contains(id)){
                d->atomIds[id] = atom;
                break;
            }
        }
    }

    return id;
}

void MoleculeEditor::setId(Atom *atom, int id)
{
    d->atomIds[id] = atom;
}

void MoleculeEditor::canUndoChangedSlot(bool canUndo)
{
    emit canUndoChanged(canUndo);
}

void MoleculeEditor::canRedoChangedSlot(bool canRedo)
{
    emit canRedoChanged(canRedo);
}

} // end chemkit namespace
