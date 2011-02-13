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

#include "displaysettingsdock.h"
#include "ui_displaysettingsdock.h"

#include "builderwindow.h"

#include <chemkit/molecule.h>
#include <chemkit/moleculewatcher.h>
#include <chemkit/graphicsmoleculeitem.h>

DisplaySettingsDock::DisplaySettingsDock(BuilderWindow *builder)
    : QDockWidget(builder),
      ui(new Ui::DisplaySettingsDock),
      m_builder(builder),
      m_showHydrogens(false)
{
    ui->setupUi(this);

    connect(ui->moleculeTypeComboBox, SIGNAL(currentIndexChanged(int)), SLOT(moleculeDisplayTypeChanged(int)));
    connect(ui->showHydrogensCheckBox, SIGNAL(clicked(bool)), SLOT(showHydrogensCheckClicked(bool)));
    connect(ui->showBondOrderCheckBox, SIGNAL(clicked(bool)), SLOT(showBondOrderCheckClicked(bool)));
    connect(builder, SIGNAL(moleculeChanged(chemkit::Molecule*)), SLOT(moleculeChanged(chemkit::Molecule*)));
}

DisplaySettingsDock::~DisplaySettingsDock()
{
    delete ui;
}

void DisplaySettingsDock::setShowHydrogens(bool showHydrogens)
{
    m_showHydrogens = showHydrogens;

    chemkit::GraphicsMoleculeItem *moleculeItem = m_builder->moleculeItem();
    if(moleculeItem)
        moleculeItem->setHydrogensVisible(showHydrogens);

    m_builder->view()->update();
}

void DisplaySettingsDock::moleculeDisplayTypeChanged(int index)
{
    chemkit::GraphicsMoleculeItem *moleculeItem = m_builder->moleculeItem();
    if(!moleculeItem)
        return;

    if(index == 0)
        moleculeItem->setDisplayType(chemkit::GraphicsMoleculeItem::BallAndStick);
    else if(index == 1)
        moleculeItem->setDisplayType(chemkit::GraphicsMoleculeItem::Stick);
    else if(index == 2)
        moleculeItem->setDisplayType(chemkit::GraphicsMoleculeItem::SpaceFilling);
}

void DisplaySettingsDock::showHydrogensCheckClicked(bool checked)
{
    setShowHydrogens(checked);
}

void DisplaySettingsDock::showBondOrderCheckClicked(bool checked)
{
    chemkit::GraphicsMoleculeItem *moleculeItem = m_builder->moleculeItem();
    if(moleculeItem)
        moleculeItem->setBondOrderVisible(checked);
}

void DisplaySettingsDock::moleculeChanged(chemkit::Molecule *molecule)
{
    Q_UNUSED(molecule);

    moleculeDisplayTypeChanged(ui->moleculeTypeComboBox->currentIndex());
}
