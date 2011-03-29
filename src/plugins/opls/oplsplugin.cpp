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

#include "oplsplugin.h"

#include "oplsatomtyper.h"
#include "oplsforcefield.h"

OplsPlugin::OplsPlugin()
    : chemkit::Plugin("opls")
{
    registerPluginClass<chemkit::AtomTyper>("opls", createOplsAtomTyper);

    chemkit::ForceField::registerForceField("opls", createOplsForceField);
}

OplsPlugin::~OplsPlugin()
{
    unregisterPluginClass<chemkit::AtomTyper>("opls");

    chemkit::ForceField::unregisterForceField("opls", createOplsForceField);
}

chemkit::AtomTyper* OplsPlugin::createOplsAtomTyper()
{
    return new OplsAtomTyper;
}

chemkit::ForceField* OplsPlugin::createOplsForceField()
{
    return new OplsForceField;
}

Q_EXPORT_PLUGIN2(opls, OplsPlugin);
