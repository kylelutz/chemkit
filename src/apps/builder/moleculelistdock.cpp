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

#include "moleculelistdock.h"
#include "ui_moleculelistdock.h"

#include "builderwindow.h"
#include "moleculepropertiesdialog.h"

MoleculeListDock::MoleculeListDock(BuilderWindow *builder)
    : QDockWidget(builder),
      ui(new Ui::MoleculeListDock)
{
    m_builder = builder;

    ui->setupUi(this);
    connect(ui->tableWidget, SIGNAL(itemSelectionChanged()), SLOT(itemSelectionChanged()));
    connect(ui->tableWidget, SIGNAL(itemDoubleClicked(QTableWidgetItem*)), SLOT(itemDoubleClicked(QTableWidgetItem*)));
    connect(ui->tableWidget, SIGNAL(itemChanged(QTableWidgetItem*)), SLOT(itemChanged(QTableWidgetItem*)));
    connect(ui->tableWidget, SIGNAL(customContextMenuRequested(const QPoint&)), SLOT(customContextMenuRequested(const QPoint&)));

    connect(builder, SIGNAL(fileChanged(chemkit::MoleculeFile*)), SLOT(fileChanged(chemkit::MoleculeFile*)));
    connect(builder, SIGNAL(moleculeChanged(chemkit::Molecule*)), SLOT(moleculeChanged(chemkit::Molecule*)));

    m_selectedItem = 0;
}

MoleculeListDock::~MoleculeListDock()
{
    delete ui;
}

void MoleculeListDock::fileChanged(chemkit::MoleculeFile *file)
{
    ui->tableWidget->clearContents();

    if(file){
        ui->tableWidget->setRowCount(file->moleculeCount());

        for(size_t i = 0; i < file->moleculeCount(); i++){
            chemkit::Molecule *molecule = file->molecule(i).get();
            QTableWidgetItem *item = new QTableWidgetItem(molecule->name().c_str());
            ui->tableWidget->setItem(i, 0, item);
        }
        ui->tableWidget->setCurrentCell(0, 0);
    }
}

void MoleculeListDock::moleculeChanged(chemkit::Molecule *molecule)
{
    Q_UNUSED(molecule);
}

void MoleculeListDock::itemSelectionChanged()
{
    boost::shared_ptr<chemkit::Molecule> molecule = currentMolecule();
    if(molecule){
        m_builder->setMolecule(molecule);
    }
}

void MoleculeListDock::itemDoubleClicked(QTableWidgetItem *item)
{
    Q_UNUSED(item);

    showMoleculeProperties();
}

void MoleculeListDock::itemChanged(QTableWidgetItem *item)
{
    boost::shared_ptr<chemkit::Molecule> molecule = currentMolecule();
    if(!molecule)
        return;

    QByteArray text = item->text().toAscii();

    if(molecule->name() != text.constData()){
        molecule->setName(text.constData());
    }
}

void MoleculeListDock::customContextMenuRequested(const QPoint &pos)
{
    QTableWidget *tableWidget = ui->tableWidget;

    QTableWidgetItem *item = tableWidget->itemAt(pos);
    if(!item)
        return;

    QMenu menu;
    QAction *renameAction = menu.addAction("Rename");
    connect(renameAction, SIGNAL(triggered()), this, SLOT(renameMolecule()));
    QAction *deleteAction = menu.addAction("Delete");
    connect(deleteAction, SIGNAL(triggered()), this, SLOT(deleteMolecule()));
    QAction *propertiesAction = menu.addAction("Properties");
    connect(propertiesAction, SIGNAL(triggered()), this, SLOT(showMoleculeProperties()));

    menu.exec(tableWidget->viewport()->mapToGlobal(pos));
}

void MoleculeListDock::renameMolecule()
{
    ui->tableWidget->editItem(ui->tableWidget->currentItem());
}

void MoleculeListDock::deleteMolecule()
{
    boost::shared_ptr<chemkit::Molecule> molecule = currentMolecule();
    if(molecule){
        m_builder->file()->removeMolecule(molecule);
    }
}

void MoleculeListDock::showMoleculeProperties()
{
    const boost::shared_ptr<chemkit::Molecule> &molecule = currentMolecule();
    if(molecule){
        MoleculePropertiesDialog dialog(molecule.get(), m_builder);
        dialog.exec();
    }
}

boost::shared_ptr<chemkit::Molecule> MoleculeListDock::currentMolecule() const
{
    if(m_builder->file()){
        int row = ui->tableWidget->currentRow();

        if(row == -1){
            return boost::shared_ptr<chemkit::Molecule>();
        }
        else if(row < static_cast<int>(m_builder->file()->moleculeCount())){
            return m_builder->file()->molecule(row);
        }
    }

    return boost::shared_ptr<chemkit::Molecule>();
}
