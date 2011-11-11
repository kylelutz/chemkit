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

#include "inchiplugin.h"

#ifdef CHEMKIT_WITH_IO
#include <chemkit/moleculefileformatadaptor.h>
#endif

#include "inchilineformat.h"
#include "inchikeylineformat.h"

InchiPlugin::InchiPlugin()
    : chemkit::Plugin("inchi")
{
    registerPluginClass<chemkit::LineFormat>("inchi", createInchiFormat);
    registerPluginClass<chemkit::LineFormat>("inchikey", createInchiKeyFormat);

#ifdef CHEMKIT_WITH_IO
    registerPluginClass<chemkit::io::MoleculeFileFormat>("inchi", createInchiFileFormat);
#endif
}

InchiPlugin::~InchiPlugin()
{
    unregisterPluginClass<chemkit::LineFormat>("inchi");
    unregisterPluginClass<chemkit::LineFormat>("inchikey");

#ifdef CHEMKIT_WITH_IO
    unregisterPluginClass<chemkit::io::MoleculeFileFormat>("inchi");
#endif
}

chemkit::LineFormat* InchiPlugin::createInchiFormat()
{
    return new InchiLineFormat;
}

chemkit::LineFormat* InchiPlugin::createInchiKeyFormat()
{
    return new InchiKeyLineFormat;
}

#ifdef CHEMKIT_WITH_IO
chemkit::io::MoleculeFileFormat* InchiPlugin::createInchiFileFormat()
{
    return new chemkit::io::MoleculeFileFormatAdaptor<chemkit::LineFormat>(new InchiLineFormat);
}
#endif

CHEMKIT_EXPORT_PLUGIN(inchi, InchiPlugin)
