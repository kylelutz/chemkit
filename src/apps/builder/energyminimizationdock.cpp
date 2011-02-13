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

#include "energyminimizationdock.h"
#include "ui_energyminimizationdock.h"

#include "builderwindow.h"
#include "energyminimizer.h"

// --- Construction and Destruction ---------------------------------------- //
EnergyMinimizationDock::EnergyMinimizationDock(BuilderWindow *builder)
    : QDockWidget(builder),
      ui(new Ui::EnergyMinimizationDock)
{
    ui->setupUi(this);
    m_builder = builder;

    connect(ui->startButton, SIGNAL(clicked()), SLOT(buttonClicked()));
    connect(builder->energyMinimizer(), SIGNAL(stateChanged(int)), SLOT(stateChanged(int)));
    connect(ui->forceFieldSelector, SIGNAL(currentIndexChanged(int)), SLOT(forceFieldChanged(int)));

    stateChanged(EnergyMinimizer::Stopped);
}

EnergyMinimizationDock::~EnergyMinimizationDock()
{
    delete ui;
}

// --- Slots --------------------------------------------------------------- //
void EnergyMinimizationDock::buttonClicked()
{
    EnergyMinimizer *optimizer = m_builder->energyMinimizer();

    if(optimizer->state() == EnergyMinimizer::Stopped)
        optimizer->start();
    else
        optimizer->stop();
}

void EnergyMinimizationDock::stateChanged(int state)
{
    EnergyMinimizer *optimizer = m_builder->energyMinimizer();

    ui->statusText->setText(optimizer->stateString());

    if(state == EnergyMinimizer::Stopped){
        ui->startButton->setText("Start");
        ui->startButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
    }
    else{
        ui->startButton->setText("Stop");
        ui->startButton->setIcon(style()->standardIcon(QStyle::SP_MediaStop));
    }

    if(state == EnergyMinimizer::UpdateReady){
        ui->energyText->setText(QString("%1 kcal/mol").arg(optimizer->energy(), 0, 'h', 1));
    }
}

void EnergyMinimizationDock::forceFieldChanged(int index)
{
    Q_UNUSED(index);

    QString forceFieldName = ui->forceFieldSelector->currentText().toLower();
    m_builder->energyMinimizer()->setForceField(forceFieldName);
}
