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
    boost::shared_ptr<chemkit::Molecule> currentMolecule() const;

private:
    Ui::MoleculeListDock *ui;
    BuilderWindow *m_builder;
    QTableWidgetItem *m_selectedItem;
};

#endif // MOLECULELISTDOCK_H
