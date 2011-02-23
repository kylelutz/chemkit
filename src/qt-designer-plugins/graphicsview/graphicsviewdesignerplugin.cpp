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

#include "graphicsviewdesignerplugin.h"

#include <chemkit/graphicsview.h>

// --- Construction and Destruction ---------------------------------------- //
GraphicsViewDesignerPlugin::GraphicsViewDesignerPlugin(QObject *parent)
    : QObject(parent)
{
}

GraphicsViewDesignerPlugin::~GraphicsViewDesignerPlugin()
{
}

// --- Properties ---------------------------------------------------------- //
QString GraphicsViewDesignerPlugin::name() const
{
    return "chemkit::GraphicsView";
}

QString GraphicsViewDesignerPlugin::group() const
{
    return "Display Widgets [chemkit]";
}

QString GraphicsViewDesignerPlugin::includeFile() const
{
    return "<chemkit/graphicsview.h>";
}

QIcon GraphicsViewDesignerPlugin::icon() const
{
    return QIcon();
}

QString GraphicsViewDesignerPlugin::toolTip() const
{
    return QString();
}

QString GraphicsViewDesignerPlugin::whatsThis() const
{
    return QString();
}

bool GraphicsViewDesignerPlugin::isContainer() const
{
    return false;
}

QWidget* GraphicsViewDesignerPlugin::createWidget(QWidget *parent)
{
    return new chemkit::GraphicsView(parent);
}

QString GraphicsViewDesignerPlugin::domXml() const
{
    return "<ui language=\"c++\">"
               "<widget class=\"chemkit::GraphicsView\" name=\"graphicsView\">"
                   "<property name=\"geometry\">"
                       "<rect>"
                           "<x>0</x>"
                           "<y>0</y>"
                           "<width>200</width>"
                           "<height>200</height>"
                       "</rect>"
                   "</property>"
               "</widget>"
           "</ui>";
}

Q_EXPORT_PLUGIN2(graphicsviewdesignerplugin, GraphicsViewDesignerPlugin)
