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

#ifndef CHEMKIT_MOLECULEEDITOR_H
#define CHEMKIT_MOLECULEEDITOR_H

#include "widgets.h"

#include <QObject>

#include <chemkit/molecule.h>

namespace chemkit {

class MoleculeEditorPrivate;

class CHEMKIT_WIDGETS_EXPORT MoleculeEditor : public QObject
{
    Q_OBJECT

public:
    // construction and destruction
    MoleculeEditor(Molecule *molecule = 0);
    ~MoleculeEditor();

    // properties
    void setMolecule(Molecule *molecule);
    Molecule* molecule() const;

    // editing
    void undo();
    bool canUndo() const;
    void redo();
    bool canRedo() const;
    void clearUndoStack();
    void beginEdit();
    void endEdit();
    bool isInEdit() const;
    void cut(const std::vector<Atom *> &atoms);
    void copy(const std::vector<Atom *> &atoms);
    std::vector<Atom *> paste();
    bool canPaste() const;
    std::vector<Atom *> copyBuffer() const;
    void clearCopyBuffer();

    // modification
    Atom* addAtom(const Element &element);
    Atom* addAtomCopy(const Atom *atom);
    void removeAtom(Atom *atom);
    void setAtomElement(Atom *atom, const Element &element);
    void setAtomPosition(Atom *atom, const Point3 &position);
    Bond* addBond(Atom *a, Atom *b, int order = 1);
    void removeBond(Bond *bond);
    void setBondOrder(Bond *bond, int order);

    // internal methods
    Atom* atom(int id);
    Bond* bond(int id1, int id2);
    int id(Atom *atom);
    void setId(Atom *atom, int id);

signals:
    void canUndoChanged(bool canUndo);
    void canRedoChanged(bool canRedo);
    void canPasteChanged(bool canPaste);

private slots:
    void canUndoChangedSlot(bool canUndo);
    void canRedoChangedSlot(bool canRedo);

private:
    MoleculeEditorPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEEDITOR_H
