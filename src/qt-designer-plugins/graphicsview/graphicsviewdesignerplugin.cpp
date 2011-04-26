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
