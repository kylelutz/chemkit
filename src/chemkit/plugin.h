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

#ifndef CHEMKIT_PLUGIN_H
#define CHEMKIT_PLUGIN_H

#include "chemkit.h"

#include <QtCore>

namespace chemkit {

class PluginPrivate;

class CHEMKIT_EXPORT Plugin : public QObject
{
    Q_OBJECT

    public:
        // properties
        QString name() const;
        QString dataPath() const;

    protected:
        // construction and destruction
        Plugin(const QString &name);
        virtual ~Plugin();

    private:
        void setFileName(const QString &fileName);

        friend class PluginManager;

    private:
        PluginPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PLUGIN_H
