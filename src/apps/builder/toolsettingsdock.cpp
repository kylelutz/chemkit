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

#include "toolsettingsdock.h"
#include "ui_toolsettingsdock.h"

#include "buildertool.h"

ToolSettingsDock::ToolSettingsDock(BuilderWindow *parent)
    : QDockWidget(parent),
      ui(new Ui::ToolSettingsDock)
{
    ui->setupUi(this);

    m_parent = parent;
    m_settingsWidget = 0;
    connect(ui->comboBox, SIGNAL(activated(int)), this, SLOT(toolComboBoxChanged(int)));
    connect(builder(), SIGNAL(toolChanged(BuilderTool *)), this, SLOT(toolChanged(BuilderTool *)));
}

ToolSettingsDock::~ToolSettingsDock()
{
    delete ui;
}

void ToolSettingsDock::toolChanged(BuilderTool *tool)
{
    if(tool == builder()->navigateTool())
        ui->comboBox->setCurrentIndex(0);
    else if(tool == builder()->buildTool())
        ui->comboBox->setCurrentIndex(1);
    else if(tool == builder()->manipulateTool())
        ui->comboBox->setCurrentIndex(2);

    if(m_settingsWidget){
        ui->layout->removeWidget(m_settingsWidget);
        delete m_settingsWidget;
    }

    m_settingsWidget = tool->settingsWidget();
    if(m_settingsWidget){
        ui->layout->addWidget(m_settingsWidget);
    }
}

void ToolSettingsDock::toolComboBoxChanged(int index)
{
    // 0 = navigate, 1 = build, 2 = manipulate
    if(index == 0)
        builder()->setTool(builder()->navigateTool());
    else if(index == 1)
        builder()->setTool(builder()->buildTool());
    else if(index == 2)
        builder()->setTool(builder()->manipulateTool());
}
