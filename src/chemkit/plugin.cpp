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

#include "plugin.h"

namespace chemkit {

// === PluginPrivate ======================================================= //
class PluginPrivate
{
    public:
        QString name;
        QString fileName;
};

// === Plugin ============================================================== //
/// \class Plugin plugin.h chemkit/plugin.h
/// \ingroup chemkit
/// \brief The Plugin class is the base class for all dynamically
///        loaded chemkit plugins.
///
/// \sa PluginManager

// --- Construction and Destruction ---------------------------------------- //
Plugin::Plugin(const QString &name)
    : QObject(),
      d(new PluginPrivate)
{
    d->name = name;
}

Plugin::~Plugin()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the plugin.
QString Plugin::name() const
{
    return d->name;
}

QString Plugin::dataPath() const
{
    return QFileInfo(d->fileName).path() + "/data/" + d->name + "/";
}

// --- Internal Methods ---------------------------------------------------- //
void Plugin::setFileName(const QString &fileName)
{
    d->fileName = fileName;
}

} // end chemkit namespace
