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

#include <QtCore>
#include <QtDesigner>

class GraphicsViewDesignerPlugin : public QObject, public QDesignerCustomWidgetInterface
{
    Q_OBJECT
    Q_INTERFACES(QDesignerCustomWidgetInterface)

    public:
        // construction and destruction
        GraphicsViewDesignerPlugin(QObject *parent = 0);
        ~GraphicsViewDesignerPlugin();

        // properties
        QString name() const;
        QString group() const;
        QString includeFile() const;
        QIcon icon() const;
        QString toolTip() const;
        QString whatsThis() const;
        bool isContainer() const;
        QWidget* createWidget(QWidget *parent);
        QString domXml() const;
};
