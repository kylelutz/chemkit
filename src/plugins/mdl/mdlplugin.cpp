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

#include "mdlplugin.h"

#include "mdlfileformat.h"

MdlPlugin::MdlPlugin()
    : chemkit::Plugin("mdl")
{
    registerPluginClass<chemkit::ChemicalFileFormat>("mdl", createMdlFormat);
    registerPluginClass<chemkit::ChemicalFileFormat>("mol", createMolFormat);
    registerPluginClass<chemkit::ChemicalFileFormat>("sdf", createSdfFormat);
    registerPluginClass<chemkit::ChemicalFileFormat>("sd", createSdFormat);
}

MdlPlugin::~MdlPlugin()
{
    unregisterPluginClass<chemkit::ChemicalFileFormat>("mdl");
    unregisterPluginClass<chemkit::ChemicalFileFormat>("mol");
    unregisterPluginClass<chemkit::ChemicalFileFormat>("sdf");
    unregisterPluginClass<chemkit::ChemicalFileFormat>("sd");
}

chemkit::ChemicalFileFormat* MdlPlugin::createMdlFormat()
{
    return new MdlFileFormat("mdl");
}

chemkit::ChemicalFileFormat* MdlPlugin::createMolFormat()
{
    return new MdlFileFormat("mol");
}

chemkit::ChemicalFileFormat* MdlPlugin::createSdfFormat()
{
    return new MdlFileFormat("sdf");
}

chemkit::ChemicalFileFormat* MdlPlugin::createSdFormat()
{
    return new MdlFileFormat("sd");
}

Q_EXPORT_PLUGIN2(mdl, MdlPlugin);
