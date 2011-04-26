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
