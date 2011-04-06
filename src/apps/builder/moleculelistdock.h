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

#ifndef MOLECULELISTDOCK_H
#define MOLECULELISTDOCK_H

#include <QtGui>

#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

namespace Ui {
    class MoleculeListDock;
}

class BuilderWindow;

class MoleculeListDock : public QDockWidget
{
    Q_OBJECT

    public:
        MoleculeListDock(BuilderWindow *builder);
        ~MoleculeListDock();

    private slots:
        void fileChanged(chemkit::MoleculeFile *file);
        void moleculeChanged(chemkit::Molecule *molecule);
        void itemSelectionChanged();
        void itemDoubleClicked(QTableWidgetItem *item);
        void itemChanged(QTableWidgetItem *item);
        void customContextMenuRequested(const QPoint &pos);
        void renameMolecule();
        void deleteMolecule();
        void showMoleculeProperties();

    private:
        chemkit::Molecule* currentMolecule() const;

    private:
        Ui::MoleculeListDock *ui;
        BuilderWindow *m_builder;
        QTableWidgetItem *m_selectedItem;
};

#endif // MOLECULELISTDOCK_H
