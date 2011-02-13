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

#ifndef TOOLSETTINGSDOCK_H
#define TOOLSETTINGSDOCK_H

#include <QDockWidget>

#include "builderwindow.h"

namespace Ui
{
    class ToolSettingsDock;
}

class ToolSettingsDock : public QDockWidget
{
    Q_OBJECT

    public:
        ToolSettingsDock(BuilderWindow *parent = 0);
        ~ToolSettingsDock();

        BuilderWindow* builder() const { return m_parent; }

    private slots:
        void toolChanged(BuilderTool *tool);
        void toolComboBoxChanged(int index);

    private:
        Ui::ToolSettingsDock *ui;
        BuilderWindow *m_parent;
        QWidget *m_settingsWidget;
};

#endif // TOOLSETTINGSDOCK_H
