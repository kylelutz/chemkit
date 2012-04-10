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

// These files implement the MMFF force field.
//
// Some useful references:
//     Description of MMFF in Towhee:
//          http://towhee.sourceforge.net/forcefields/mmff94.html
//     Parameter Description from CHARMM:
//          http://www.charmm.org/documentation/c32b2/mmff_params.html
//     MMFF Validation Suite:
//          http://server.ccl.net/cca/data/MMFF94/
//     Parameter Data Files:
//          ftp://ftp.wiley.com/public/journals/jcc/suppmat/17/490/MMFF-I_AppendixB.ascii

#include "mmffforcefield.h"

#include "mmffatomtyper.h"
#include "mmffparameters.h"
#include "mmffcalculation.h"

#include <chemkit/plugin.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/topology.h>
#include <chemkit/pluginmanager.h>

// --- Construction and Destruction ---------------------------------------- //
MmffForceField::MmffForceField()
    : chemkit::ForceField("mmff"),
      m_parameters(0)
{
    const chemkit::Plugin *mmffPlugin = chemkit::PluginManager::instance()->plugin("mmff");
    if(mmffPlugin){
        std::string dataPath = mmffPlugin->dataPath();
        addParameterSet("mmff94", dataPath + "mmff94.prm");
        setParameterSet("mmff94");
    }

    setFlags(chemkit::ForceField::AnalyticalGradient);
}

MmffForceField::~MmffForceField()
{
    delete m_parameters;
}

// --- Parameterization ---------------------------------------------------- //
bool MmffForceField::setup()
{
    if(!m_parameters || m_parameters->fileName() != parameterFile()){
        m_parameters = new MmffParameters;
        bool ok = m_parameters->read(parameterFile());
        if(!ok){
            setErrorString("Failed to load parameters: " + m_parameters->errorString());
            delete m_parameters;
            m_parameters = 0;
            return false;
        }
    }

    const boost::shared_ptr<chemkit::Topology> &topology = this->topology();
    if(!topology){
        return false;
    }

    // bond strech calculations
    foreach(const chemkit::Topology::BondedInteraction &interaction, topology->bondedInteractions()){
        size_t a = interaction[0];
        size_t b = interaction[1];

        addCalculation(new MmffBondStrechCalculation(a, b));
    }

    // angle bend and strech bend calculations
    foreach(const chemkit::Topology::AngleInteraction &interaction, topology->angleInteractions()){
        size_t a = interaction[0];
        size_t b = interaction[1];
        size_t c = interaction[2];

        addCalculation(new MmffAngleBendCalculation(a, b, c));
        addCalculation(new MmffStrechBendCalculation(a, b, c));
    }

    // out of plane bending calculation (for each trigonal center)
    foreach(const chemkit::Topology::ImproperTorsionInteraction &interaction, topology->improperTorsionInteractions()){
        size_t a = interaction[0];
        size_t b = interaction[1];
        size_t c = interaction[2];
        size_t d = interaction[3];

        addCalculation(new MmffOutOfPlaneBendingCalculation(a, b, c, d));
        addCalculation(new MmffOutOfPlaneBendingCalculation(a, b, d, c));
        addCalculation(new MmffOutOfPlaneBendingCalculation(c, b, d, a));
    }

    // torsion calculations (for each dihedral)
    foreach(const chemkit::Topology::TorsionInteraction &interaction, topology->torsionInteractions()){
        size_t a = interaction[0];
        size_t b = interaction[1];
        size_t c = interaction[2];
        size_t d = interaction[3];

        addCalculation(new MmffTorsionCalculation(a, b, c, d));
    }

    // van der waals and electrostatic calculations
    foreach(const chemkit::Topology::NonbondedInteraction &interaction, topology->nonbondedInteractions()){
        size_t a = interaction[0];
        size_t b = interaction[1];

        addCalculation(new MmffVanDerWaalsCalculation(a, b));
        addCalculation(new MmffElectrostaticCalculation(a, b));
    }

    bool ok = true;

    foreach(chemkit::ForceFieldCalculation *calculation, calculations()){
        bool setup = static_cast<MmffCalculation *>(calculation)->setup(m_parameters);

        if(!setup){
            ok = false;
        }

        setCalculationSetup(calculation, setup);
    }

    return ok;
}

const MmffParameters* MmffForceField::parameters() const
{
    return m_parameters;
}
