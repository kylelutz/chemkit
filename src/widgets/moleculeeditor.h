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

#ifndef CHEMKIT_MOLECULEEDITOR_H
#define CHEMKIT_MOLECULEEDITOR_H

#include "widgets.h"

#include <QObject>
#include <QUndoStack>

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
        void cut(const QList<Atom *> &atoms);
        void copy(const QList<Atom *> &atoms);
        QList<Atom *> paste();
        bool canPaste() const;
        QList<Atom *> copyBuffer() const;
        void clearCopyBuffer();

        // modification
        Atom* addAtom(const Element &element);
        Atom* addAtomCopy(const Atom *atom);
        void removeAtom(Atom *atom);
        void setAtomAtomicNumber(Atom *atom, int atomicNumber);
        void setAtomPosition(Atom *atom, const Point3 &position);
        Bond* addBond(Atom *a, Atom *b, int order = Bond::Single);
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
