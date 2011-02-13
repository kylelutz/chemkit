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

    connect(builder, SIGNAL(fileChanged(chemkit::ChemicalFile*)), SLOT(fileChanged(chemkit::ChemicalFile*)));
    connect(builder, SIGNAL(moleculeChanged(chemkit::Molecule*)), SLOT(moleculeChanged(chemkit::Molecule*)));

    m_selectedItem = 0;
}

MoleculeListDock::~MoleculeListDock()
{
    delete ui;
}

void MoleculeListDock::fileChanged(chemkit::ChemicalFile *file)
{
    ui->tableWidget->clearContents();

    if(file){
        ui->tableWidget->setRowCount(file->moleculeCount());

        for(int i = 0; i < file->moleculeCount(); i++){
            chemkit::Molecule *molecule = file->molecule(i);
            QTableWidgetItem *item = new QTableWidgetItem(molecule->name());
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
    chemkit::Molecule *molecule = currentMolecule();
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
    chemkit::Molecule *molecule = currentMolecule();
    if(!molecule)
        return;

    if(molecule->name() != item->text()){
        molecule->setName(item->text());
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
    chemkit::Molecule *molecule = currentMolecule();
    if(molecule){
        m_builder->file()->removeMolecule(molecule);
    }
}

void MoleculeListDock::showMoleculeProperties()
{
    const chemkit::Molecule *molecule = currentMolecule();
    if(molecule){
        MoleculePropertiesDialog dialog(molecule, m_builder);
        dialog.exec();
    }
}

chemkit::Molecule* MoleculeListDock::currentMolecule() const
{
    if(m_builder->file()){
        int row = ui->tableWidget->currentRow();
        if(row < m_builder->file()->moleculeCount()){
            return m_builder->file()->molecule(ui->tableWidget->currentRow());
        }
    }

    return 0;
}
