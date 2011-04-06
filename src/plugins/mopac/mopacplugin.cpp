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

#include "mopacplugin.h"

#include "mopinfileformat.h"
#include "mopcrtfileformat.h"

MopacPlugin::MopacPlugin()
    : chemkit::Plugin("mopac")
{
    registerPluginClass<chemkit::MoleculeFileFormat>("mopin", createMopinFormat);
    registerPluginClass<chemkit::MoleculeFileFormat>("mopcrt", createMopcrtFormat);
}

MopacPlugin::~MopacPlugin()
{
    unregisterPluginClass<chemkit::MoleculeFileFormat>("mopin");
    unregisterPluginClass<chemkit::MoleculeFileFormat>("mopcrt");
}

chemkit::MoleculeFileFormat* MopacPlugin::createMopinFormat()
{
    return new MopinFileFormat;
}

chemkit::MoleculeFileFormat* MopacPlugin::createMopcrtFormat()
{
    return new MopcrtFileFormat;
}

Q_EXPORT_PLUGIN2(mopac, MopacPlugin);
