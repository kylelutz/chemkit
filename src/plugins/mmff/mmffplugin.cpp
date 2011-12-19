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

#include "mmffplugin.h"

#include <chemkit/forcefieldenergydescriptor.h>

#include "mmffatomtyper.h"
#include "mmffforcefield.h"
#include "mmffparametersdata.h"
#include "mmffaromaticitymodel.h"
#include "mmffpartialchargepredictor.h"

MmffPlugin::MmffPlugin()
    : chemkit::Plugin("mmff")
{
    registerPluginClass<chemkit::AtomTyper>("mmff", createMmffAtomTyper);
    registerPluginClass<chemkit::ForceField>("mmff", createMmffForceField);
    registerPluginClass<chemkit::AromaticityModel>("mmff", createMmffAromaticityModel);
    registerPluginClass<chemkit::MolecularDescriptor>("mmff-energy", createMmffEnergyDescriptor);
    registerPluginClass<chemkit::PartialChargePredictor>("mmff", createMmffPartialChargePredictor);
}

MmffPlugin::~MmffPlugin()
{
    foreach(MmffParametersData *parameters, m_parametersCache.values()){
        parameters->deref();
    }

    unregisterPluginClass<chemkit::AtomTyper>("mmff");
    unregisterPluginClass<chemkit::ForceField>("mmff");
    unregisterPluginClass<chemkit::AromaticityModel>("mmff");
    unregisterPluginClass<chemkit::PartialChargePredictor>("mmff");
}

void MmffPlugin::storeParameters(const QString &name, MmffParametersData *parameters)
{
    if(m_parametersCache.contains(name)){
        m_parametersCache[name]->deref();
    }

    m_parametersCache.insert(name, parameters);
    parameters->ref();
}

MmffParametersData* MmffPlugin::parameters(const QString &name) const
{
    return m_parametersCache.value(name, 0);
}

chemkit::AtomTyper* MmffPlugin::createMmffAtomTyper()
{
    return new MmffAtomTyper;
}

chemkit::ForceField* MmffPlugin::createMmffForceField()
{
    return new MmffForceField;
}

chemkit::AromaticityModel* MmffPlugin::createMmffAromaticityModel()
{
    return new MmffAromaticityModel;
}

chemkit::MolecularDescriptor* MmffPlugin::createMmffEnergyDescriptor()
{
    return new chemkit::ForceFieldEnergyDescriptor<MmffForceField>("mmff-energy");
}

chemkit::PartialChargePredictor* MmffPlugin::createMmffPartialChargePredictor()
{
    return new MmffPartialChargePredictor;
}

CHEMKIT_EXPORT_PLUGIN(mmff, MmffPlugin)
