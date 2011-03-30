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

#include "inchiplugin.h"

#include "inchilineformat.h"
#include "inchikeylineformat.h"
#include "inchifileformat.h"

InchiPlugin::InchiPlugin()
    : chemkit::Plugin("inchi")
{
    registerPluginClass<chemkit::LineFormat>("inchi", createInchiFormat);
    registerPluginClass<chemkit::LineFormat>("inchikey", createInchiKeyFormat);
    registerPluginClass<chemkit::ChemicalFileFormat>("inchi", createInchiFileFormat);
}

InchiPlugin::~InchiPlugin()
{
    unregisterPluginClass<chemkit::LineFormat>("inchi");
    unregisterPluginClass<chemkit::LineFormat>("inchikey");
    unregisterPluginClass<chemkit::ChemicalFileFormat>("inchi");
}

chemkit::LineFormat* InchiPlugin::createInchiFormat()
{
    return new InchiLineFormat;
}

chemkit::LineFormat* InchiPlugin::createInchiKeyFormat()
{
    return new InchiKeyLineFormat;
}

chemkit::ChemicalFileFormat* InchiPlugin::createInchiFileFormat()
{
    return new InchiFileFormat;
}

Q_EXPORT_PLUGIN2(inchi, InchiPlugin);
